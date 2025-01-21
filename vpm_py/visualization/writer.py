from matplotlib.animation import FFMpegWriter
import gc
from contextlib import contextmanager

class OptimizedFFMpegWriter(FFMpegWriter):
    '''FFMpeg-pipe writer bypassing figure.savefig with memory optimizations.'''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.frame_format = 'argb'
        self._cached_w = None
        self._cached_h = None
        self._cached_dpi = None

    def setup(self, fig, outfile, dpi=None):
        super().setup(fig, outfile, dpi)
        # Cache initial dimensions
        self._cached_w = self._w
        self._cached_h = self._h
        self._cached_dpi = self.dpi
        self.frame_count = 0
        
    @contextmanager
    def _frame_context(self):
        try:
            yield
        finally:
            # Clear figure and force garbage collection
            gc.collect()

    def grab_frame(self, **savefig_kwargs):
        '''Memory-optimized frame grabber'''
        with self._frame_context():
            try:
                # Only adjust size if changed
                if (self.fig.get_figwidth() != self._cached_w or 
                    self.fig.get_figheight() != self._cached_h):
                    
                    self.fig.set_size_inches(self._cached_w, self._cached_h)

                if (self.fig.dpi != self._cached_dpi):
                    self.fig.set_dpi(self._cached_dpi)
                
                # Write frame and immediately flush
                self._proc.stdin.write(self.fig.canvas.tostring_argb())
                self._proc.stdin.flush()
                self.frame_count += 1
                
            except (RuntimeError, IOError) as e:
                out, err = self._proc.communicate()
                raise IOError('Error saving animation to file (cause: {0}) '
                            'Stdout: {1} StdError: {2}. It may help to re-run '
                            'with --verbose-debug.'.format(e, out, err))
    
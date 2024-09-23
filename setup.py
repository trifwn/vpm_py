# -*- coding: utf-8 -*-
#   Copyright (c) 2019 D. de Vries
#
import os
import platform
import re
import subprocess
import sys

from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext

__version__ = re.findall(
    r"""__version__ = ["']+([0-9\.]*)["']+""",
    open('vpm_py/__init__.py').read(),
)[0]

options = {k: 'OFF' for k in ['--opt', '--debug', '--cuda']}
for flag in options.keys():
    if flag in sys.argv:
        options[flag] = 'ON'
        sys.argv.remove(flag)

# Command line flags forwarded to CMake
cmake_cmd_args = []
for f in sys.argv:
    if f.startswith('-D'):
        cmake_cmd_args.append(f)
        sys.argv.remove(f)


class CMakeExtension(Extension):

    def __init__(self, name, cmake_list_dir='.', **kwargs):
        super().__init__(name, sources=[], **kwargs)
        self.cmake_lists_dir = os.path.abspath(cmake_list_dir)


class CMakeBuild(build_ext):
    def run(self):
        # Ensure that CMake is present and working
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError('Cannot find CMake executable')
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        extdir = os.path.join(extdir, "shared_libs")
        os.makedirs(extdir, exist_ok=True)
        cfg ='{}'.format('Release' if self.debug else 'Debug')

        cmake_args = [
            f'-DCMAKE_BUILD_TYPE={cfg}',
            # Ask CMake to place the resulting library in the directory
            # containing the extension
            f'-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}',
            # Other intermediate static libraries are placed in a
            # temporary build directory instead
            f'-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{cfg.upper()}={self.build_temp}',
            # Hint CMake to use the same Python executable that is launching the build, prevents possible
            # mismatching if multiple versions of Python are installed
            '-DPYTHON_EXECUTABLE={}'.format(sys.executable),
        ]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''), self.distribution.get_version())
        
        os.makedirs(self.build_temp, exist_ok=True)

        # Config and build the extension
        # Print the command in a way that can be copied and pasted
        print(f'Configuring extension for {ext.name}. Command Run is:\n')
        print(' '.join(['cmake', ext.cmake_lists_dir] + cmake_args))
        print('-'* 10)

        subprocess.check_call(
            ['cmake', ext.cmake_lists_dir] + cmake_args,
            cwd=self.build_temp, env=env
        )
        
        print(f'Building extension for {ext.name}. Command Run is:\n')
        print(' '.join(['cmake', '--build', '.', '--config', cfg]))
        print('-'* 10)
        build_args = ['--config', f'{cfg}']
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)


def readme():
    return "Hello"

setup(
    name='vpm_py',
    version=__version__,
    description='compiled python module ',
    long_description=readme(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering',
    ],
    keywords='vpm aerodynamic fluid mechanics',
    author='Tryfonas Themas',
    author_email='tryfonthem@gmail.com',
    packages=['vpm_py'],
    ext_modules=[CMakeExtension('vpm_py.vpm')],
    cmdclass={'build_ext': CMakeBuild},
    install_requires=[
        # Numerics
        'numpy', 
        "scipy",

        # Visualization
        "matplotlib",
        "pyQt5",

        # Data
        "h5py",
        "pandas",
        "tqdm",

        # Parallelism
        "mpi4py", # CONDA-FORGE

        # MATH
        "mkl", # CONDA
        # "mkl-service", # CONDA

        # Intel libraries
        'intel-fortran-rt', 
        'intel-cmplr-lib-rt',
        'intel-cmplr-lic-rt',
        "dpcpp-cpp-rt",
        ],
    zip_safe=False
)

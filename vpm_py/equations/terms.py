
class CoordinateTerm:
    def __init__(self, name: str, n_components: int):
        self.name = name
        self.n_components = n_components

    def __str__(self):
        return self.name

    def to_latex(self):
        return r"\text{" + self.name + "}"
    
    def get_dimensionality(self):
        return (self.n_components,)

    def get_structure(self):
        return {self.name : self}



from . import Expression, DifferentialExpression, CoordinateTerm

class DifferentialOperator(Expression):
    def __init__(self, name: str, latex_repr: str, input_dim: int, output_dim: int):
        super().__init__()
        self.name = name
        self.latex_repr = latex_repr
        self.input_dim = input_dim
        self.output_dim = output_dim

    def __call__(self, term: Expression, basis: CoordinateTerm):
        if self.input_dim is None:
            self.input_dim = term.get_dimensionality()
        if self.output_dim is None:
            self.output_dim = term.get_dimensionality()
        return DifferentialExpression(self, [term, basis])
    
    def __str__(self):
        return self.name

    def to_latex(self):
        return self.latex_repr

    def get_dimensionality(self):
        return self.output_dim

    def get_structure(self):
        return {self.name : self}

class Gradient(DifferentialOperator):
    def __init__(self):
        super().__init__("∇", r"\nabla", input_dim=None, output_dim=None)

    def __call__(self, term: Expression, basis: CoordinateTerm):
        if term.get_dimensionality() != 0:
            raise ValueError("Gradient operator can only be applied to scalar fields.")
        scalar_field = term
        self.input_dim = 1
        self.output_dim = basis.n_components
        return super().__call__(scalar_field, basis)

class Divergence(DifferentialOperator):
    def __init__(self):
        super().__init__("∇⋅", r"\nabla \cdot", input_dim=None, output_dim=None)

    def __call__(self, term: Expression, basis: CoordinateTerm):
        if term.get_dimensionality() != basis.n_components:
            raise ValueError("Dimension mismatch between vector field and basis.")
        vector_field = term        
        self.input_dim = vector_field.get_dimensionality()
        self.output_dim = 1  # Divergence is a scalar
        return super().__call__(vector_field, basis)

class Curl(DifferentialOperator):
    def __init__(self):
        super().__init__("∇×", r"\nabla \times", input_dim=None, output_dim=None)

    def __call__(self, term: Expression, basis: CoordinateTerm): 
        if term.get_dimensionality() != basis.n_components:
            raise ValueError("Dimension mismatch between vector field and basis.")
        vector_field = term
        self.input_dim = vector_field.get_dimensionality()
        self.output_dim = vector_field.get_dimensionality()  # Curl output dimensionality
        return super().__call__(vector_field, basis)

class Laplacian(DifferentialOperator):
    def __init__(self):
        super().__init__("∇²", r"\nabla^2", input_dim=None, output_dim=None)

    def __call__(self, term: Expression, basis: CoordinateTerm):
        if term.get_dimensionality() != 0:
            raise ValueError("Laplacian operator can only be applied to scalar fields.")
        scalar_field = term
        self.input_dim = 1
        self.output_dim = 1
        return super().__call__(scalar_field, basis)

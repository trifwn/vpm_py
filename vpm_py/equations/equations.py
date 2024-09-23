from abc import ABC, abstractmethod
from . import Expression, AlgerbraicOperator, Constant

class Equation:
    num_equations = 0

    def __init__(self, lhs: Expression | float | int, rhs: Expression | float | int, name: str | None = None):
        self.lhs: Expression = Constant(str(lhs)) if isinstance(lhs, (int, float)) else lhs
        self.rhs: Expression = Constant(str(rhs)) if isinstance(rhs, (int, float)) else rhs
        self._validate_dimensions()
        Equation.num_equations += 1
        self.name = name if name else f"Equation {Equation.num_equations}"

    def _validate_dimensions(self):
        lhs_dim = self.lhs.get_dimensionality()
        rhs_dim = self.rhs.get_dimensionality()
        if lhs_dim != rhs_dim and lhs_dim != 0 and rhs_dim != 0:
            raise ValueError("Dimensionality mismatch between LHS and RHS in the equation.")

    def __str__(self):
        return f"{self.name}: {self.lhs} = {self.rhs}"

    def to_latex(self):
        return f"{self.lhs.to_latex()} = {self.rhs.to_latex()}"

    def extract_unknowns(self):
        unknowns = []
        if isinstance(self.rhs, AlgerbraicOperator):
            unknowns.extend(self.rhs.get_terms())
        return unknowns

    def get_structure(self):
        return {
            self.name : self 
        }

class Problem:
    def __init__(self, name: str, equations: list[Equation]):
        self.name = name
        self.equations = equations

    def __str__(self):
        return f"{self.name}:\n" + "\n".join(str(eq) for eq in self.equations)

    def to_latex(self):
        return r"\text{" + self.name + r"}:\\" + r"\\".join(eq.to_latex() for eq in self.equations)

    def extract_unknowns(self):
        unknowns = set()
        for eq in self.equations:
            unknowns.update(eq.extract_unknowns())
        return list(unknowns)

    def get_problem_structure(self):
        structure = []
        for eq in self.equations:
            structure.append(eq.get_structure())
        return {self : structure}

    def is_well_formulated(self):
        return len(self.extract_unknowns()) == len(self.equations)

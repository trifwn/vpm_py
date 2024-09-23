from __future__ import annotations
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from . import DifferentialOperator, CoordinateTerm

class Expression(ABC):
    def __init__(self):
        self.terms = []  
        self.operations = []

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def to_latex(self):
        pass

    @abstractmethod
    def get_dimensionality(self):
        pass

    def operation(self, op: str, other: 'Expression'):
        # Check dimensionality compatibility
        left_dim = self.get_dimensionality()
        right_dim = other.get_dimensionality()

        if left_dim != right_dim and left_dim != 0 and right_dim != 0:
            raise ValueError(f"Dimension mismatch: {left_dim} and {right_dim}")

        if left_dim == 0:
            left_dim = right_dim

        # Store operation relationship
        self.terms.append((self, op, other))
        self.operations.append(op)

        return AlgerbraicOperator(self, op, other, result_dim=left_dim)

    def get_terms(self):
        return [term for term in self.terms]

    def get_structure(self):
        structure = []
        for term in self.terms:
            if isinstance(term, AlgerbraicOperator):
                structure.append(term.get_structure())
            else:
                structure.append(str(term))
        return structure

    def __add__(self, other):
        return self.operation("+", other)

    def __sub__(self, other):
        return self.operation("-", other)

    def __mul__(self, other):
        return self.operation("*", other)

    def __truediv__(self, other):
        return self.operation("/", other)

    def __neg__(self):
        return AlgerbraicOperator(Constant("-1"), "*", self, result_dim=self.get_dimensionality())

    def __pow__(self, other):
        return self.operation("^", other)

class AlgerbraicOperator(Expression):
    def __init__(self, left: Expression, op: str, right: Expression, result_dim: int):
        super().__init__()
        self.left = left
        self.op = op
        self.right = right
        self.result_dim = result_dim
        self.terms.append((left, op, right))

    def get_dimensionality(self):
        return self.result_dim

    def __str__(self):
        return f"({self.left} {self.op} {self.right})"

    def to_latex(self):
        op_map = {"+": "+", "-": "-", "*": r"\cdot", "/": r"\frac{#1}{#2}", "^": "^"}
        if self.op == "/":
            return op_map[self.op].replace("#1", f"{{{self.left.to_latex()}}}").replace("#2", f"{{{self.right.to_latex()}}}")
        return f"({self.left.to_latex()} {op_map[self.op]} {self.right.to_latex()})"

    def get_structure(self):
        return {self.op : self}

class PDETerm(Expression):
    _instances: dict[tuple[str, tuple[str], int], 'PDETerm'] = {}

    def __new__(cls, name: str, function_of: list, num_dimensions: int = 0):
        key = (name, tuple(function_of), num_dimensions)
        if key in cls._instances:
            return cls._instances[key]
        instance = super().__new__(cls)
        cls._instances[key] = instance
        return instance

    def __init__(self, name: str, function_of: list, num_dimensions: int):
        super().__init__()
        self.name = name
        self.function_of = function_of
        self.num_dimensions = num_dimensions

    def __str__(self):
        return self.name

    def to_latex(self):
        return r"\text{" + self.name + "}"

    def get_dimensionality(self):
        return self.num_dimensions
    
    def get_structure(self):
        return {self.name : self}

class Constant(PDETerm):
    _constant_instances: dict[str, 'Constant'] = {}

    def __new__(cls, name: str):
        key = name
        if key in cls._constant_instances:
            return cls._constant_instances[key]
        instance = super().__new__(cls, name=name, function_of=[], num_dimensions=0)
        cls._constant_instances[key] = instance
        return instance

    def __init__(self, name: str):
        super().__init__(name=name, function_of=[], num_dimensions=0)

    def __str__(self):
        return str(self.name)

    def to_latex(self):
        return str(self.name)

    def get_dimensionality(self):
        return 0

    def get_structure(self):
        return {self.name : self}

class DifferentialExpression(Expression):
    def __init__(self, operator: DifferentialOperator, terms: list[Expression | CoordinateTerm]):
        super().__init__()
        self.operator = operator
        self.terms = terms

    def __str__(self):
        return f"{self.operator}({', '.join(map(str, self.terms))})"

    def to_latex(self):
        terms_latex = ", ".join(term.to_latex() for term in self.terms)
        return f"{self.operator.to_latex()}({terms_latex})"

    def get_dimensionality(self):
        return self.operator.get_dimensionality()

    def get_structure(self):
        return {self.operator : self}

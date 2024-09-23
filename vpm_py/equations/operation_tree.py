from . import (
    Problem, 
    Equation, 
    Constant, 
    PDETerm, 
    CoordinateTerm, 
    DifferentialOperator, 
    DifferentialExpression, 
    AlgerbraicOperator
)

class TreeNode:
    def __init__(self, value):
        self.value = value
        self.num_leaf_nodes = 0
        self.leaf_nodes = []

    def __repr__(self):
        return f"TreeNode({self.value})"
    
    def __str__(self):
        return f"TreeNode({self.value}), {self.num_leaf_nodes} leaf nodes"

class OperationTree:
    def __init__(self):
        self.root = None

    def build_tree(self, structure):
        self.root = self._build_tree_recursive(structure)

    def _build_tree_recursive(self, node):
        if node is None:
            return None
        return self.decide_based_on_key(node)

    def decide_based_on_key(self, key):
        if isinstance(key, (str, int, float)):
            return TreeNode(str(key))
        elif isinstance(key, dict):
            return self.add_dict_node(key)
        elif isinstance(key, (list, tuple)):
            return self.add_list_or_tuple_node(key)
        elif isinstance(key, Problem):
            return self.add_problem_node(key)
        elif isinstance(key, Equation):
            return self.add_equation_node(key)
        elif isinstance(key, Constant):
            return self.add_constant_node(key)
        elif isinstance(key, PDETerm):
            return self.add_pde_term_node(key)
        elif isinstance(key, CoordinateTerm):
            return self.add_coordinate_term_node(key)
        elif isinstance(key, DifferentialOperator):
            return self.add_differential_operator_node(key)
        elif isinstance(key, DifferentialExpression):
            return self.add_differential_expression_node(key)
        elif isinstance(key, AlgerbraicOperator):
            return self.add_operation_node(key)
        else:
            raise ValueError(f"Unknown type: {type(key)}")

    def add_dict_node(self, key):
        nodes = []
        for k, v in key.items():
            nodes.append(self.decide_based_on_key(k))
        return nodes[0] if len(nodes) == 1 else nodes

    def add_list_or_tuple_node(self, key):
        nodes = []
        for item in key:
            nodes.append(self.decide_based_on_key(item))
        return nodes[0] if len(nodes) == 1 else nodes

    def add_problem_node(self, key):
        node = TreeNode(f"Problem: {key.name}")
        for eq in key.equations:
            node.leaf_nodes.append(self.decide_based_on_key(eq))
        node.num_leaf_nodes = len(node.leaf_nodes)
        return node

    def add_equation_node(self, key):
        node = TreeNode(f"Equation: {key.name}\n=")
        node.leaf_nodes.append(self.decide_based_on_key(key.lhs))
        node.leaf_nodes.append(self.decide_based_on_key(key.rhs))
        node.num_leaf_nodes = len(node.leaf_nodes)
        return node

    def add_pde_term_node(self, key):
        node = TreeNode(f"{key.name}({', '.join(map(str, key.function_of))})")
        return node

    def add_constant_node(self, key):
        return TreeNode(f"{key.name}")

    def add_coordinate_term_node(self, key):
        node = TreeNode(f"{key.name}")
        return node

    def add_differential_operator_node(self, key):
        node = TreeNode(f"{key.name}")
        return node

    def add_differential_expression_node(self, key: DifferentialExpression):
        operator = key.operator
        node = TreeNode(f"{operator}")
        terms = key.terms
        for term in terms:
            node.leaf_nodes.append(self.decide_based_on_key(term))
        node.num_leaf_nodes = len(node.leaf_nodes)
        return node

    def add_operation_node(self, key):
        node = TreeNode(f"{key.op}")
        node.leaf_nodes.append(self.decide_based_on_key(key.left))
        node.leaf_nodes.append(self.decide_based_on_key(key.right))
        node.num_leaf_nodes = len(node.leaf_nodes)
        return node

    def visualize_str(self):
        self._visualize_recursive_str(self.root, "", True)

    def _visualize_recursive_str(self, node, prefix, is_last):
        if node is not None:
            try:
                print(prefix + ("└── " if is_last else "├── ") + str(node.value))
                for i, child in enumerate(node.leaf_nodes):
                    new_prefix = prefix + ("    " if is_last else "│   ")
                    self._visualize_recursive(child, new_prefix, i == len(node.leaf_nodes) - 1)
            except AttributeError:
                print(prefix + ("└── " if is_last else "├── ") + str(node))
                return
    
    def visualize(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.axis('off')
        
        if self.root is not None:
            self._visualize_recursive_plot(self.root, ax, pos=(0, 0), level=0, width=2**(self._get_tree_height(self.root)))

        fig.tight_layout()
        fig.show()

    def _visualize_recursive_plot(self, node, ax, pos, level, width):
        if node is not None:
            # Draw the node
            ax.text(pos[0], pos[1], str(node.value), ha='center', va='center', bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='lightgray'))

            # Calculate positions for children
            child_pos = self._get_child_positions(pos, level, width)

            # Draw edges and recursively draw children
            for i, child in enumerate(node.leaf_nodes):
                if child is not None:
                    ax.plot([pos[0], child_pos[i][0]], [pos[1] - 0.2, child_pos[i][1] + 0.2], color='black')  # Draw edge
                    self._visualize_recursive_plot(child, ax, child_pos[i], level + 1, width // 2)

    def _get_child_positions(self, parent_pos, level, width):
        # Calculate positions for children based on the parent's position
        positions = []
        num_children = len(self.root.leaf_nodes)
        for i in range(num_children):
            x_offset = (i - (num_children - 1) / 2) * (width / num_children)
            positions.append((parent_pos[0] + x_offset, parent_pos[1] - 1))
        return positions

    def _get_tree_height(self, node):
        if node is None:
            return 0
        return 1 + max(self._get_tree_height(child) for child in node.leaf_nodes) if node.leaf_nodes else 0

    def flatten(self):
        return self._flatten_recursive(self.root)
    
    def _flatten_recursive(self, node):
        if node is not None:
            try:
                return [node.value] + [self._flatten_recursive(child) for child in node.leaf_nodes]
            except AttributeError:
                return [node]
        return []

from mpi4py import MPI
import numpy as np
from vpm_py.operators_lib import OperatorsLib

# Initialize MPI and OperatorsLib
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
operators = OperatorsLib()

# Define grid parameters
dx = 1e-3
dy = 1e-3
dz = 1e-3
x = np.arange(-50*dx, 50*dx, dx)
y = np.arange(-50*dy, 50*dy, dy)
z = np.arange(-50*dz, 50*dz, dz)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

def print_green(text):
    print(f"\033[92m{text}\033[00m")

def print_red(text):
    print(f"\033[91m{text}\033[00m")

test_number = 1
def run_test(operator_name, field, analytical_result, calc_function, *args):
    global test_number
    if rank == 0:
        calc_start_time = MPI.Wtime()
        result_np = calc_function(field, dx, dy, dz, *args)
        calc_end_time = MPI.Wtime()

        print(f"\t{test_number}:\t{operator_name} Test Results:")
        test_number += 1
        print(f"\t\tCalc Time: {calc_end_time - calc_start_time:.6f} seconds")
        error = np.abs(result_np - analytical_result)
        print(f"\t\tMAX Error: {error.max():.6e}")
        print(f"\t\tError (mean): {error.mean():.6e}")
        print(f"\t\tError (std): {error.std():.6e}")
        if np.max(error) < 1e-5:
            print_green("\tAll tests passed.")
        else:
            print_red("\tSome tests failed.")
        print()

#############################################################################
# Scalar Field: f = sin(X) * cos(Y) * Z^2
#############################################################################
scalar_field =  np.sin(X) * np.cos(Y) * Z**2

scalar_df_dx = np.cos(X) * np.cos(Y) * Z**2
scalar_df_dy = -np.sin(X) * np.sin(Y) * Z**2
scalar_df_dz = 2 * np.sin(X) * np.cos(Y) * Z

scalar_df_dx2 = -np.sin(X) * np.cos(Y) * Z**2     # ∂²f/∂x²
scalar_df_dy2 = -np.sin(X) * np.cos(Y) * Z**2     # ∂²f/∂y²
scalar_df_dz2 = 2.0 * np.sin(X) * np.cos(Y)       # ∂²f/∂z²

scalar_df_dxdy= -np.cos(X) * np.sin(Y) * Z**2     # ∂²f/∂x∂y
scalar_df_dxdz= 2.0 * np.cos(X) * np.cos(Y) * Z   # ∂²f/∂x∂z
scalar_df_dydz= -2.0 * np.sin(X) * np.sin(Y) * Z  # ∂²f/∂y∂z

#############################################################################
# Vector Field: f = [X^2 + Y^2 + Z^2, X*Y + Y*Z + Z*X, sin(X)*cos(Y)*Z^2]
#############################################################################
f1 = X**2 + Y**2 + Z**2
f2 = X*Y + Y*Z + Z*X
f3 = np.sin(X) * np.cos(Y) * Z**2
vector_field = np.stack([f1, f2, f3])

df1_dx = 2 * X
df1_dy = 2 * Y
df1_dz = 2 * Z

df2_dx = Y + Z
df2_dy = X + Z  
df2_dz = X + Y

df3_dx = np.cos(X) * np.cos(Y) * Z**2
df3_dy = -np.sin(X) * np.sin(Y) * Z**2
df3_dz = 2 * np.sin(X) * np.cos(Y) * Z  

# f3 = sin(x) * cos(y) * z^2
df3_dx2 = -np.sin(X) * np.cos(Y) * Z**2     # ∂²f3/∂x²
df3_dy2 = -np.sin(X) * np.cos(Y) * Z**2     # ∂²f3/∂y²
df3_dz2 = 2.0 * np.sin(X) * np.cos(Y)       # ∂²f3/∂z²
df3_dxdy= -np.cos(X) * np.sin(Y) * Z**2     # ∂²f3/∂x∂y
df3_dxdz= 2.0 * np.cos(X) * np.cos(Y) * Z   # ∂²f3/∂x∂z
df3_dydz= -2.0 * np.sin(X) * np.sin(Y) * Z  # ∂²f3/∂y∂z

# f1 = x^2 + y^2 + z^2
df1_dx2 = 2.0  * np.ones_like(X)           # ∂²f1/∂x²
df1_dy2 = 2.0  * np.ones_like(X)           # ∂²f1/∂y²
df1_dz2 = 2.0  * np.ones_like(X)           # ∂²f1/∂z²
df1_dxdy = 0.0 * np.ones_like(X)           # ∂²f1/∂x∂y
df1_dxdz = 0.0 * np.ones_like(X)           # ∂²f1/∂x∂z
df1_dydz = 0.0 * np.ones_like(X)           # ∂²f1/∂y∂z

# f2 = x*y + y*z + z*x
df2_dx2 =  0.0 * np.ones_like(X)            # ∂²f2/∂x²
df2_dy2 =  0.0 * np.ones_like(X)            # ∂²f2/∂y²
df2_dz2 =  0.0 * np.ones_like(X)            # ∂²f2/∂z²
df2_dxdy = 1.0 * np.ones_like(X)            # ∂²f2/∂x∂y
df2_dxdz = 1.0 * np.ones_like(X)            # ∂²f2/∂x∂z
df2_dydz = 1.0 * np.ones_like(X)            # ∂²f2/∂y∂z
#############################################################################

#############################################################################
# Scalar Field
#############################################################################

# Test 1: Derivative (of a scalar field)
print("Running Derivative Test")
run_test("Derivative", scalar_field, scalar_df_dx, operators.calc_derivative, 1, 1) 

# Test 2: Gradient (of a scalar field)
gradient_field = np.stack([scalar_df_dx, scalar_df_dy, scalar_df_dz])
print("Running Gradient Test")
run_test("Gradient", scalar_field, gradient_field, operators.calc_gradient)

# Test 3: Laplacian (of a scalar field)
laplacian_field = scalar_df_dx2 + scalar_df_dy2 + scalar_df_dz2
print("Running Scalar Field Laplacian Test")
run_test("Laplacian", scalar_field, laplacian_field, operators.calc_laplacian)

# Test 4: Jacobian (of a scalar field) (same as Gradient)
jacobian_field = np.stack([scalar_df_dx, scalar_df_dy, scalar_df_dz])
print("Running Scalar Field Jacobian Test")
run_test("Jacobian", scalar_field, jacobian_field, operators.calc_jacobian)

# Test 5: Hessian (of a scalar field)
hessian_field = np.stack([
    [scalar_df_dx2, scalar_df_dxdy, scalar_df_dxdz],
    [scalar_df_dxdy, scalar_df_dy2, scalar_df_dydz],
    [scalar_df_dxdz, scalar_df_dydz, scalar_df_dz2]
]) 
print("Running Scalar Field Hessian Test")
run_test("Hessian", scalar_field, hessian_field, operators.calc_hessian)


#############################################################################
# Vector Field
#############################################################################

# Test 6: Derivative (of a vector field) with respect to X
df_field = np.stack([df1_dx, df2_dx, df3_dx])
print("Running Vector Derivative Test")
run_test("Derivative", vector_field, df_field, operators.calc_derivative, 1,1)

# Test 7: Divergence (of a vector field)
div_field = df1_dx + df2_dy + df3_dz  # Sum of the components for divergence
print("Running Divergence Test")
run_test("Divergence", vector_field, div_field, operators.calc_divergence)

# Test 8: Curl (of a vector field)
curl_field = np.stack([
    df3_dy - df2_dz,  # d(f3)/dy - d(f2)/dz
    df1_dz - df3_dx,  # d(f1)/dz - d(f3)/dx
    df2_dx - df1_dy   # d(f2)/dx - d(f1)/dy
])
print("Running Curl Test")
run_test("Curl", vector_field, curl_field, operators.calc_curl)

# Test 9: Vector Laplacian
lap_f1 = df1_dx2 + df1_dy2 + df1_dz2
lap_f2 = df2_dx2 + df2_dy2 + df2_dz2
lap_f3 = df3_dx2 + df3_dy2 + df3_dz2
vector_laplacian_field = np.stack([lap_f1, lap_f2, lap_f3])
print("Running Vector Laplacian Test")
run_test("Vector Laplacian", vector_field, vector_laplacian_field, operators.calc_vector_laplacian)

# Test 10: Jacobian matrix 
jacobian_field = np.stack([
    [df1_dx, df1_dy, df1_dz],
    [df2_dx, df2_dy, df2_dz],
    [df3_dx, df3_dy, df3_dz]
])
print("Running Jacobian Test")
run_test("Jacobian", vector_field, jacobian_field, operators.calc_jacobian)

# Test 11: Hessian (of a vector field)
hessian_f1 = np.stack([ # Hessian of f1 = X^2 + Y^2 + Z^2
    [df1_dx2,  df1_dxdy, df1_dxdz],
    [df1_dxdy, df1_dy2,  df1_dydz],
    [df1_dxdz, df1_dydz, df1_dz2 ]
])
hessian_f2 = np.stack([ # Hessian of f2 = X*Y + Y*Z + Z*X
    [df2_dx2,  df2_dxdy, df2_dxdz],
    [df2_dxdy, df2_dy2,  df2_dydz],
    [df2_dxdz, df2_dydz, df2_dz2 ]
])
hessian_f3 = np.stack([ # Hessian of f3 = sin(X)*cos(Y)*Z^2
    [df3_dx2,  df3_dxdy, df3_dxdz],
    [df3_dxdy, df3_dy2,  df3_dydz],
    [df3_dxdz, df3_dydz, df3_dz2 ]
])
hessian_field = np.stack([hessian_f1, hessian_f2, hessian_f3])
print("Running Vector Hessian Test")
run_test("Vector Hessian", vector_field, hessian_field, operators.calc_hessian)

if rank == 0:
    print("All tests completed.")
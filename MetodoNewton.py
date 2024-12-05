import sympy as sp

#Definir variables de las funciones

x = sp.symbols('x')
y = sp.symbols('y')

#Ingresar funciones a utilizar

#funciones = [x**2 + y**2 - 25,sp.exp(x) + y - 15]
funciones = [x**2 + y**2 - 30,sp.sin(x) + y - 10]

#Econtrar derivadas en las funciones

derivadas = []
for i in range (len(funciones)):
        df_dx = sp.diff(funciones[i],x)
        df_dy = sp.diff(funciones[i],y)
        derivadas.append([df_dx,df_dy])

#Resolución del sistema

def newton(funciones,valor,valor_fx,derivadas):
        iteraciones = 1
        error_relativo = sp.Matrix([10000,10000])
        valor_ant = sp.Matrix(valor)
        derivadas_val = [[0,0],[0,0]]
        
        while max(error_relativo) > 0.001:
                valor_fx.clear()
                #Valorzar funcion F(x)
                for i in range (len(funciones)):
                    resultado = funciones[i].subs({x:valor[i],y:valor[i]}).evalf()
                    valor_fx.append(resultado)

                #Valorizar funcion J(x)
                for i in range(len(derivadas)):
                    for j in range(len(derivadas)):
                        derivadas_val[i][j] =  derivadas[i][j].subs({x:valor[i],y:valor[i]}).evalf()

                derivadas_val_matrix = sp.Matrix(derivadas_val)
        
                # Verificar si el determinante es cero
                determinante = derivadas_val_matrix.det()
                if determinante == 0:
                    print("Error: La matriz Jacobiana no es invertible (determinante = 0).")
                    print(f"Valores actuales: {valor}")
                    print(f"Iteraciones realizadas: {iteraciones}")
                    return
        
                # Calcular la inversa de J(x)
                derivadas_inv = derivadas_val_matrix.inv()

                #Calcular delta x
                delta_x = -(derivadas_inv * sp.Matrix(valor_fx))

                #Calculamos el valor de la siguite iteración
                valor = sp.Matrix(valor) + delta_x

                for i in range (len(valor)):
                        error_relativo[i] = abs(delta_x[i]/valor[i])
                print(error_relativo)

                iteraciones = iteraciones + 1 
        print(f"Valores de X: {valor}")
        print(f"Numero de iteraciones: {iteraciones}")
        print(f"Error relativo: {error_relativo}")


valor = [1,1]
valor_fx = []
newton(funciones, valor, valor_fx,derivadas)
print(f"Funciones: {funciones}")
print(f"Derivadas de J(x): {derivadas}")

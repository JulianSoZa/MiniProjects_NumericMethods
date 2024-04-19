En este trabajo se desarrollan los problemas propuestos en el miniproyecto numero 2 de la asignatura de Métodos Numéricos para Ingeniería, en el cual se aborda el método de diferencias finitas para solucionar el potencial y el campo electrico en coordenadas polares y cartesianas en 2D para un capacitor en forma de anillo.

El trabajo consta de modulos con funciones para resolver un problema específico y un archivo main donde se podrá correr la rutina que contiene los procedimientos, o puede comentar alguna parte en específica.

A continuacion se describe la funcion de cada módulo, tenga en cuenta que podrá probar por separado su funcionamineto dentro de estos mismos:

El paquete domainDiscretization contiene dos modulos:
        El modulo cartesian se utiliza para hacer la discretizacion del dominio en coordenadas cartesianas.
        El modulo polar se utiliza para hacer la discretizacion en coordenadas polares.

El modulo electric_potencial tiene las funciones:
        electric_potential_solution: para calcular el potencial electrico en coordenadas polares.
        electric_potential_solution_cartesian: para calcular el potencial electrico en coordenadas cartesianas. 

El modulo electric_field tiene las funciones:
        electric_field_solution: para calcular el campo electrico en coordenadas polares
        electric_field_solution_cartesian: para calcular el campo electrico en coordenadas cartesianas

El modulo plot_electric_solution tiene las funciones:
        ploter_finite_solutions: para graficar la solucion del potencial y el campo electrico en coordenadas polares.
        ploter_finite_solutions_cartesian: para graficar la solucion del potencial y el campo electrico en coordenadas cartesianas. 

El modulo analysis_results para obtener las gráficas de la comparacion y el analisis de resultados.
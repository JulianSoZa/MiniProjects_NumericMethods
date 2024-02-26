En este cuaderno se desarrollan los problemas propuestos en el miniproyecto numero 1 de la asignatura de Métodos Numéricos para Ingeniería, en el cual se abordaron diferentes problemas de los método numéricos como la aproximación e interpolación de funciones, la integración numérica 1D y 2D, y el tratamiento de mallas.

El cuaderno consta de cinco secciones en donde se resuelven los cinco enunciados del miniproyecto.

A continuación de describe el procedimiento a seguir para cada una de las secciones:

Antes de ejecutar el codigo que contienen las secciones, ejecute la celda en la que se especifica las librerias o modulos que se utilizan. Por favor asegúrese que tiene instaladas las librerías y los modulos correspondientes. 

Sección 1. Método de proyección de Galerkin

Sección 3. Integración 1D: --------------------------------------

    Ejecute la celda en la que se ecuentran las funciones para la integración numérica:
            ----- numerical_integration_1D
            ----- segmented_1D_numerical_integration

    Luego, para tener el valor de referencia a comparar, ejecute la celda en donde se cálcula la solución analítica de la integral.
            ------ En esta puede modificar la función que desea integrar y los límites de integración. Y se puede observar el resultado de la integral analítica en consola.

    Finalmente, la ultima celda de la sección corresponde al cálculo de la integral numérica de las dos formas.
            ------ En esta puede variar el parámetro n de los dos métodos, que significa su valor final en la iteración para la realizacion de la grafica de error relativo vs el parámetro n. Cabe recordar que el parámetro n inicia con n = 2.

            ------ Si quiere conocer el valor de la integral para cada valor de n en cada uno de los métodos, puede imprimir en consola los arreglos:
                ----- n_Integrals_SegmentedDomain
                ----- n_Integrals_FullDomain

            ------ Si quiere conocer el valor del error relativo para cada valor de n en cada uno de los métodos, puede imprimir en consola los arreglos:
                ----- errorR_SegmentedDomain
                ----- errorR_FullDomain

Sección 4. Mallas e Integracion 2D: -----------------------------
    
    Para tener el valor de referencia a comparar, ejecute la celda en donde se cálcula la solución analítica de la integral.
            ------ En esta puede modificar la función que desea integrar y los límites de integración (siempre y cuando estos últimos correspondan a un dominio rectangular). Y se puede observar el resultado de la integral analítica en consola.

    Seguidamente, ejecute la celda en la que se encuentra la función **gaussian_integration_2D** creada por Daniel Giraldo Cuartas.

    Luego, ejecute la celda en donde se almacena en una malla y se grafica el dominio original dado.
            ------ En esta puede modificar la malla inicial dada siempre y cuando esté compuesta por cuadrilateros y se definan los vertices. 

    A continuación ejecute la celda en la que se procede a subdividir en rectangulos la malla anterior.
            ------ En esta puede modificar los parámetros:
                ----- "n" : numero de elementos a lo largo del eje horizontal.
                ----- "z" : numero de elementos a lo largo del eje vertical. Si quiere que los elementos sean cuadrados n debe ser multiplo de 5 y utilice la linea en la que z depende de n.
            ------ Se puede obervar la malla subdividida.

    Finalmente ejecute la celda en la que se calcula la integral numérica para cada elemento de la malla.
            ------ En esta puede modificar la función a integrar y el numero de puntos de cuadratura para los dos ejes 'x' e 'y'.
            ------ Se puede observar en consola la integral numérica para todo el domino y el error relativo con repecto a la integral anlítica.
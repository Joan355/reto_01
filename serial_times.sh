#!/bin/bash

# Nombre del ejecutable de C
ejecutable="./jacobi1d"

# Ciclo para iterar sobre la cantidad de hilos

    # Nombre del archivo CSV donde guardar los resultados
    archivo_resultados="resultados_hilos.csv"

    # Cabecera del archivo CSV
    echo "N,Nodos,Tiempo" > $archivo_resultados

    # Ciclo para iterar sobre la cantidad de nodos
    for nodos in {100..1000..200}; do
        # Ciclo para iterar sobre el rango de N
        for n in {500..10000..500}; do
            # Ejecutar el programa C con los argumentos N, nodos y hilos
            tiempo=$(./jacobi1d $n $nodos)

            # Agregar el resultado al archivo CSV
            echo "$n,$nodos,$tiempo" >> $archivo_resultados
        done
    done



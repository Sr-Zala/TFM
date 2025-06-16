#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#Script de preprocesamiento FASTQ sin dependencias externas.


import os
import sys
import shutil
import time
import gzip
from pathlib import Path
from collections import defaultdict
import json

# Adaptadores conocidos
ADAPTADORES = {
    "TruSeq": "AGATCGGAAGAG",
    "SmallRNA": "TGGAATTCTCGG",
    "Nextera": "CTGTCTCTTATA"
}

# Lectura de archivos FASTQ, línea a línea
def leer_fastq(archivo):
    opener = gzip.open if archivo.suffix == '.gz' else open
    with opener(archivo, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            secuencia = f.readline().strip()
            _ = f.readline()  # línea +
            calidades = f.readline().strip()
            yield header, secuencia, calidades

# Conversión de carácter a valor Phred
def phred_a_valor(caracter):
    return ord(caracter) - 33

# Aplicación de los filtros sobre una secuencia individual
def procesar_secuencia(header, secuencia, calidades, params):
    if params['adaptador'] != 'auto':
        secuencia = secuencia.replace(params['adaptador'], '')
    else:
        for adaptador in ADAPTADORES.values():
            secuencia = secuencia.replace(adaptador, '')
    
    secuencia = secuencia[params['trim_left']:]
    if params['trim_right'] > 0:
        secuencia = secuencia[:-params['trim_right']]
    calidades = calidades[params['trim_left']:]
    if params['trim_right'] > 0:
        calidades = calidades[:-params['trim_right']]
    
    valores_phred = [phred_a_valor(q) for q in calidades]
    if not valores_phred:
        return None
    elif min(valores_phred) < params['min_quality']:
        return None
    elif sum(valores_phred)/len(valores_phred) < params['min_mean_quality']:
        return None
    elif not (params['min_length'] <= len(secuencia) <= (params['max_length'] or float('inf'))):
        return None
    elif secuencia.count('N') > params['max_n']:
        return None
    else:
        return f"{header}\n{secuencia}\n+\n{calidades}\n"

# Procesa un archivo completo aplicando todos los filtros
def procesar_archivo(archivo_entrada, params, dir_salida):
    secuencias_vistas = set()
    total = 0
    conservadas = 0
    nombre_salida = f"procesado_{archivo_entrada.stem}.fastq"
    ruta_salida = dir_salida / nombre_salida
    
    with open(ruta_salida, 'w') as out:
        for header, seq, qual in leer_fastq(archivo_entrada):
            total += 1
            if params['remove_duplicates']:
                if seq in secuencias_vistas:
                    continue
                secuencias_vistas.add(seq)
            sec_procesada = procesar_secuencia(header, seq, qual, params)
            if sec_procesada:
                out.write(sec_procesada)
                conservadas += 1
    
    return {
        'archivo': archivo_entrada.name,
        'total': total,
        'conservadas': conservadas,
        'porcentaje': (conservadas/total)*100 if total > 0 else 0
    }

# Solicita parámetros desde JSON o por consola
def obtener_parametros():
    ruta_parametros = Path(__file__).parent / "parametros.json"
    if ruta_parametros.exists():
        with open(ruta_parametros, "r") as f:
            params = json.load(f)
            if params.get("max_length") == 0:
                params["max_length"] = None
            return params

    print("\n" + "="*50)
    print("PARÁMETROS DE PROCESAMIENTO".center(50))
    print("="*50)

    params = {
        'min_quality': int(input("\nCalidad mínima por base (Phred 0-40): ")),
        'min_mean_quality': int(input("Calidad media mínima: ")),
        'min_length': int(input("Longitud mínima (bp): ")),
        'max_length': int(input("Longitud máxima (0=sin límite): ")) or None,
        'trim_left': int(input("Bases a recortar en 5': ")),
        'trim_right': int(input("Bases a recortar en 3': ")),
        'max_n': int(input("Máximo de bases 'N' permitidas: "))
    }
    resp = input("¿Eliminar duplicados? (S/n): ").strip().lower()
    params['remove_duplicates'] = resp != 'n'
    
    if input("\n¿Detección automática de adaptadores? (s/n): ").lower() == 's':
        params['adaptador'] = 'auto'
    else:
        params['adaptador'] = input("Secuencia del adaptador: ").upper()

    return params

# Programa principal
if __name__ == "__main__":
    print("\n" + "="*60)
    print("PREPROCESAMIENTO FASTQ (Sin Biopython / FastQC)".center(60))
    print("="*60)
    
    dir_actual = Path(__file__).parent.resolve()
    dir_procesadas = dir_actual / "secuencias_procesadas"
    dir_procesadas.mkdir(exist_ok=True)
    
    archivos = [f for f in dir_actual.iterdir()
                if f.is_file() and (
                    f.name.endswith('.fastq') or
                    f.name.endswith('.fq') or
                    f.name.endswith('.fastq.gz') or
                    f.name.endswith('.fq.gz')
                )]
    
    if not archivos:
        print("\nNo se encontraron archivos FASTQ (.fastq, .fq o .gz)")
        sys.exit(1)
    
    print(f"\nArchivos encontrados ({len(archivos)}):")
    for i, archivo in enumerate(archivos, 1):
        print(f" {i}. {archivo.name}")
    
    params = obtener_parametros()
    
    resultados = []
    for archivo in archivos:
        print(f"\nProcesando {archivo.name}...")
        inicio = time.time()
        
        resultado = procesar_archivo(archivo, params, dir_procesadas)
        tiempo = time.time() - inicio
        
        resultados.append({
            **resultado,
            'tiempo': f"{tiempo:.2f}s"
        })
    
    print("\n" + "="*60)
    print("RESUMEN DE PROCESAMIENTO".center(60))
    print("="*60)
    print("\nArchivo            | Total  | Conservadas | %     | Tiempo")
    print("-"*60)
    for res in resultados:
        print(f"{res['archivo'][:15]:<17} | {res['total']:6} | {res['conservadas']:10} | {res['porcentaje']:5.1f}% | {res['tiempo']:7}")
    
    print("\n" + "="*60)
    print("Proceso completado".center(60))
    print(f"Resultados en: {dir_procesadas}".center(60))
    print("="*60)
    
    resumen_path = dir_procesadas / "resumen_procesamiento.txt"
    with open(resumen_path, "w", encoding="utf-8") as resumen_file:
        resumen_file.write("="*60 + "\n")
        resumen_file.write("RESUMEN DE PROCESAMIENTO".center(60) + "\n")
        resumen_file.write("="*60 + "\n\n")
        resumen_file.write("Archivo            | Total  | Conservadas | %     | Tiempo\n")
        resumen_file.write("-"*60 + "\n")
        for res in resultados:
            resumen_file.write(f"{res['archivo'][:15]:<17} | {res['total']:6} | {res['conservadas']:10} | {res['porcentaje']:5.1f}% | {res['tiempo']:7}\n")
        resumen_file.write("\n" + "="*60 + "\n")
        resumen_file.write("Proceso completado".center(60) + "\n")
        resumen_file.write(f"Resultados en: {dir_procesadas}".center(60) + "\n")
        resumen_file.write("="*60 + "\n")

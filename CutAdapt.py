#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#Script automatizado de preprocesamiento con Cutadapt.


import subprocess
import time
from pathlib import Path
import sys

# Parámetros definidos para aplicar los filtros
PARAMS = {
    "calidad": "20",
    "media": "25",  # No aplicable en Cutadapt
    "min_long": "50",
    "recorte_5": "5",
    "recorte_3": "5",
    "max_n": "2",    # No aplicable directamente
    "adaptador": "AGATCGGAAGAG"
}

# Mensaje de inicio
print("\n" + "="*60)
print("PROCESAMIENTO FASTQ CON CUTADAPT (Automatizado)".center(60))
print("="*60)

# Configurar rutas de entrada y salida
dir_actual = Path(__file__).parent
dir_salida = dir_actual / "secuencias_procesadas_cutadapt"
dir_salida.mkdir(exist_ok=True)

# Buscar archivos FASTQ 
archivos = [
    f for f in dir_actual.iterdir()
    if f.is_file() and f.suffix.lower() in (".fastq", ".fq", ".gz")
]
if not archivos:
    print("\nNo se encontraron archivos FASTQ (.fastq, .fq o .gz)")
    sys.exit(1)

print(f"\nArchivos encontrados: {len(archivos)}")
for i, archivo in enumerate(archivos, 1):
    print(f" {i}. {archivo.name}")

# Ejecutar Cutadapt sobre cada archivo
resultados = []
for archivo in archivos:
    salida = dir_salida / f"cutadapt_{archivo.stem}.fastq"
    print(f"\nProcesando: {archivo.name}")

    comando = [
        "cutadapt",
        "-q", PARAMS['calidad'],
        "--cut", PARAMS['recorte_5'],
        "--cut", f"-{PARAMS['recorte_3']}",
        "-m", PARAMS['min_long'],
        "-a", PARAMS['adaptador'],
        "-o", str(salida),
        str(archivo)
    ]

    t0 = time.time()
    proceso = subprocess.run(comando, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t1 = time.time()
    duracion = t1 - t0

    resultados.append({
        "archivo": archivo.name,
        "salida": salida.name,
        "tiempo": f"{duracion:.2f}s",
        "exito": proceso.returncode == 0
    })

# Mostrar resumen
print("\n" + "="*60)
print("RESUMEN DE PROCESAMIENTO CUTADAPT".center(60))
print("="*60)
print("\nArchivo original     | Archivo salida           | Tiempo   | Éxito")
print("-"*60)
for r in resultados:
    print(f"{r['archivo'][:20]:<20} | {r['salida'][:22]:<22} | {r['tiempo']:>8} | {str(r['exito'])}")

print("\n" + "="*60)
print("Procesamiento completado con Cutadapt".center(60))
print(f"Resultados en: {dir_salida}".center(60))
print("="*60)

# Guardar resumen en archivo
resumen_path = dir_salida / "resumen_cutadapt.txt"
with open(resumen_path, "w", encoding="utf-8") as f:
    f.write("="*60 + "\n")
    f.write("RESUMEN DE PROCESAMIENTO CUTADAPT".center(60) + "\n")
    f.write("="*60 + "\n\n")
    f.write("Archivo original     | Archivo salida           | Tiempo   | Éxito\n")
    f.write("-"*60 + "\n")
    for r in resultados:
        f.write(f"{r['archivo'][:20]:<20} | {r['salida'][:22]:<22} | {r['tiempo']:>8} | {str(r['exito'])}\n")
    f.write("\n" + "="*60 + "\n")
    f.write("Procesamiento completado con Cutadapt".center(60) + "\n")
    f.write(f"Resultados en: {dir_salida}".center(60) + "\n")
    f.write("="*60 + "\n")

print(f"\nResumen guardado en: {resumen_path}")

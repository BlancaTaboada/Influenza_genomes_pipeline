# Pipeline de Genomas de Influenza (limpieza → mejor referencia S4/S6 → consenso)

Este repositorio contiene un flujo en **dos etapas** para procesar lecturas Illumina y generar **consensos por segmento** para Influenza A/B. La etapa 2 elige automáticamente la **mejor referencia** utilizando el soporte de mapeo en **segmento 4 (HA)** y **segmento 6 (NA)**, y luego llama consenso y variantes.

## Scripts incluidos
- `scripts/preProcess_list_v2.sh` — Limpieza por lotes con `fastp` (apto para SGE).
- `scripts/ConsensoInfluenza_v3.sh` — Mapeo multi-referencia **mínimo**, selección de la mejor referencia por **S4/S6 (HA/NA)** y generación de consenso con `ivar`.

---

## Requisitos
- Linux (probado en clusters SGE).
- Herramientas (ajusta a tu entorno): `fastp`, `bowtie2`, `samtools`, `bcftools`, `bamtools`, `ivar`.
- Variables SGE utilizadas: `$SGE_TASK_ID`, `$NSLOTS`.

### Estructura de entradas
- FASTQ/FASTQ.GZ pareados en un directorio, con convención `R1_` / `R2_` (ej. `SampleX_R1_001.fastq.gz`, `SampleX_R2_001.fastq.gz`).
- Un archivo de lista `FileNames.txt` con **R1 en la línea N** y **R2 en la línea N+1** (así se emparejan con `$SGE_TASK_ID`). Hay un ejemplo en `examples/FileNames.txt`.

### Referencias
Coloca en `refs/` todas las referencias a evaluar e indexa cada una con Bowtie2:
```bash
bowtie2-build A_H1N1.fasta A_H1N1
bowtie2-build A_H3N2.fasta A_H3N2
bowtie2-build B_Victoria.fasta B_Victoria
bowtie2-build B_Yamagata.fasta B_Yamagata
```
Para cada referencia incluye:
- `{ref}.fasta` + su índice (`{ref}.fasta.fai`).
- `{ref}.bed` con regiones/segmentos (usar nombres consistentes para **S4_HA** y **S6_NA**).
- `{ref}.gff3` por segmento si tu flujo los requiere.

> La **selección de mejor referencia** se basa en el número de alineamientos sobre S4/S6 en cada mapeo. Asegúrate que los nombres de segmento en FASTA/BED/GFF3 coincidan con lo que el script espera (por defecto `S4_HA` y `S6_NA`).

---

## Etapa 1 — Limpieza de lecturas
Ejemplo (ajusta rangos/colas/rutas según tu cluster):
```bash
qsub -R y -l h_rt=23:59:59 -pe thread 4 -t 1-100   scripts/preProcess_list_v2.sh   /ruta/a/FASTQ   FileNames.txt   /ruta/a/Results
```
**Salida:** pares `{Sample}_Q.fastq.gz` en `Results` y logs con conteos pre/post.

## Etapa 2 — Consenso + variantes (mejor referencia por S4/S6)
Ejecuta apuntando al directorio con los FASTQ limpios y tu `FileNames.txt`:
```bash
qsub -R y -l h_rt=23:59:59 -pe thread 4 -t 1-100   scripts/ConsensoInfluenza_v3.sh   /ruta/a/Results   FileNames.txt   /ruta/absoluta/a/refs
```
**Resumen por muestra:**
1. Mapea contra **cada** referencia (`bowtie2 --very-sensitive-local`).
2. **Cuenta alineamientos** en `S4_HA` y `S6_NA` y **elige la referencia** con mayor soporte.
3. Genera BAM ordenado/indexado, VCF (`bcftools`) y consensos + variantes por segmento (`ivar consensus` / `ivar variants`).
4. Organiza salidas en `FASTA/`, `SAM/`, `BAM/`, `VIRUS/` y crea un FASTA concatenado `{sample}_All_ivar2.fasta`.

---

## Consejos y solución de problemas
- Si no se elige referencia, revisa coincidencia exacta de nombres de segmentos en FASTA/BED/GFF3.
- Para lecturas muy cortas (<60 nt), podrías ajustar `--score-min` de Bowtie2.
- `ivar consensus` por defecto usa `-q 20`, `-m 5`, `-t 0` (enmascara baja cobertura con `N`). Ajusta según tu criterio.
- Verifica que `$NSLOTS` esté definido por tu cola; si no, ajusta hilos en los scripts.

---

## Estructura del repositorio
```
influenza-genomes-pipeline/
├─ scripts/
│  ├─ preProcess_list_v2.sh
│  └─ ConsensoInfluenza_v3.sh
├─ refs/
├─ examples/
│  └─ FileNames.txt
├─ docs/
│  └─ PROTOCOL.md
└─ README.md
```

## Publicar en GitHub
```bash
git init
git add .
git commit -m "Initial commit: influenza genomes pipeline (v3 minimal)"
git branch -M main
git remote add origin https://github.com/<TU_USUARIO>/influenza-genomes-pipeline.git
git push -u origin main
```

## 🧬 Acknowledgements
Este trabajo ha sido apoyado por la Universidad Nacional Autónoma de México (UNAM) mediante el proyecto **PAPIIT-DGAPA-IN230523** otorgado
 a **Blanca Taboada**, y por la **Secretaría de Educación, Ciencia, Tecnología e Innovación de la Ciudad de México (SECTEI)** mediante el
 proyecto **SECTEI/138/2024** otorgado a **Selene Zárate**.
 
 "# Influenza_genomes_pipeline" 
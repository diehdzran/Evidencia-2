# 1. Cargar librerías
library(Biostrings)
library(ape)
library(kmer)

# 2. Definir archivos y etiquetas (Asegúrate de que la ruta sea correcta)
files <- c(
  "data/PZ357334.1.fasta", 
  "data/NC_045512.2.fasta", 
  "data/PZ055857.1.fasta",
  "data/PV125869.1.fasta", 
  "data/PZ253834.1.fasta", 
  "data/PX471304.1.fasta",
  "data/PX404732.1.fasta", 
  "data/LC929404.1.fasta", 
  "data/PV910446.1.fasta",
  "data/PZ011000.1.fasta"
)

# Nombres exactos para las etiquetas del árbol
paises <- c("USA", "China", "India", "Francia", "Alemania", 
            "Brasil", "Corea", "Japon", "Italia", "UK")

# 3. Leer y convertir las secuencias
dna_sequences <- readDNAStringSet(files)

# Convertimos el objeto DNAStringSet a una lista de vectores de caracteres
sequences_list <- as.character(dna_sequences)
sequences_list <- strsplit(sequences_list, "")
names(sequences_list) <- paises

# 3.1 Calcular longitud de las secuencias
sequence_lengths <- width(dna_sequences)
names(sequence_lengths) <- paises

# 3.2 Gráfica comparativa de número de bases
barplot(sequence_lengths, 
        main = "Comparativa del Tamaño Genómico de Variantes SARS-CoV-2",
        xlab = "País / Variante", 
        ylab = "Número de Bases (pb)", 
        col = "skyblue", 
        las = 2, 
        ylim = c(29000, 30000), 
        xpd = FALSE) 

# 4. Cálculo de k-mers (Tetranucleótidos, k=4)
counts <- kcount(sequences_list, k = 4)

# 5. Generar la matriz de distancia y el cluster
dist_matrix <- dist(counts)
hc <- hclust(dist_matrix, method = "complete")


# 6. Graficar
plot(as.dendrogram(hc), 
     main = "Árbol Filogenético Variantes de SARS-CoV-2 a nivel mundial",
     ylab = "Distancia Genética",
     xlab = "", # Dejamos vacío aquí para usar mtext abajo
     sub = "")

# Añadimos los títulos del eje X
mtext("País", side = 1, line = 3, font = 1)
mtext("Basado en frecuencias de tetranucleótidos (k-mers)", side = 1, line = 4, cex = 0.8)

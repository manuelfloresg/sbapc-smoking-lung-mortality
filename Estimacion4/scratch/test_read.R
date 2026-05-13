cat(">>> Probando lectura de _runtime_setup.R...\n")
txt <- readLines("runs/_runtime_setup.R", n = 5)
cat(">>> Lectura exitosa. Primeras líneas:\n")
print(txt)

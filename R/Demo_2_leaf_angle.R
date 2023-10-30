# Evaluate the leaf angle probability density function (PDF) from the beta distribution
# using the mean and standard deviation of the leaf inclination angle.
# The leaf angle PDF is calculated for 9 angle classes between 5 and 85 degrees in increments of 10 degrees.

# Define the leaf type
leaf <- "Planophile"
# leaf <- "Erectophile"
# leaf <- "Plagiophile"
# leaf <- "Uniform"
# leaf <- "Spherical"

# Parameters for the leaf angle distribution
lad_ave <- 26.76 * (pi/180)
lad_std <- 18.5068 * (pi/180)

# Convert to p and q parameters for the beta distribution
num <- 1 - (lad_std^2 + lad_ave^2) / (lad_ave * pi / 2)
den <- (lad_std^2 + lad_ave^2) / (lad_ave^2) - 1
p <- num / den
q <- ((pi/2) / lad_ave - 1) * p

# Calculate leaf inclination angle probability density function (PDF) and fractional abundance (lad) for 9 10-degree bins
dangle <- 10 * (pi/180)  # Leaf inclination angle increment (radians)
angle <- c(5, 15, 25, 35, 45, 55, 65, 75, 85)  # Leaf inclination angle (degrees)
angle <- angle * (pi/180)  # degrees -> radians

# Initialize vectors for beta distribution PDF and fractional abundance
beta_pdf <- numeric(length(angle))
beta_lad <- numeric(length(angle))

# Loop through each angle
for (i in 1:length(angle)) {
  x <- angle[i] / (pi/2)
  fp <- x^(p - 1)
  fq <- (1 - x)^(q - 1)
  beta_pdf[i] <- 2 / pi * fp * fq / beta(p, q)  # Leaf angle probability density function
  beta_lad[i] <- beta_pdf[i] * dangle  # Fraction of leaves in this angle bin
}

# Calculate the known solution
exact_pdf <- numeric(length(angle))
exact_lad <- numeric(length(angle))

# Loop through each angle
for (i in 1:length(angle)) {
  switch(leaf,
         'Planophile' = {
           exact_pdf[i] <- 2 / pi * (1 + cos(2 * angle[i]))
         },
         'Erectophile' = {
           exact_pdf[i] <- 2 / pi * (1 - cos(2 * angle[i]))
         },
         'Plagiophile' = {
           exact_pdf[i] <- 2 / pi * (1 - cos(4 * angle[i]))
         },
         'Uniform' = {
           exact_pdf[i] <- 2 / pi
         },
         'Spherical' = {
           exact_pdf[i] <- sin(angle[i])
         })
  # Exact relative leaf angle distribution (fraction)
  exact_lad[i] <- exact_pdf[i] * dangle
}

# Print out fractional abundance and compare with the known solution
beta_sum <- sum(beta_lad)
beta_ave <- sum(angle * beta_lad)
exact_sum <- sum(exact_lad)
exact_ave <- sum(angle * exact_lad)

cat("\n")
cat("Leaf type = ", leaf, "\n")
cat("     Angle         beta            exact \n")
for (i in 1:length(angle)) {
  cat(formatC(angle[i] * 180 / pi, digits = 6), " ", formatC(beta_lad[i], digits = 4), " ", formatC(exact_lad[i], digits = 4), "\n")
}

cat("\n")
cat("beta distribution \n")
cat("Sum of leaf angle distribution = ", formatC(beta_sum, digits = 4), "\n")
cat("Mean leaf angle = ", formatC(beta_ave * 180 / pi, digits = 4), "\n")

cat("\n")
cat("Exact solution \n")
cat("Sum of leaf angle distribution = ", formatC(exact_sum, digits = 4), "\n")
cat("Mean leaf angle = ", formatC(exact_ave * 180 / pi, digits = 4), "\n")

# Calculate Ross index
F1 <- sum(beta_lad[1:3])
F2 <- sum(beta_lad[4:6])
F3 <- sum(beta_lad[7:9])
beta_xl <- 0.5 * (abs(0.134 - F1) + abs(0.366 - F2) + abs(0.5 - F3))
if ((0.5 - F3) < 0) {
  beta_xl <- -beta_xl
}

F1 <- sum(exact_lad[1:3])
F2 <- sum(exact_lad[4:6])
F3 <- sum(exact_lad[7:9])
exact_xl <- 0.5 * (abs(0.134 - F1) + abs(0.366 - F2) + abs(0.5 - F3))
if ((0.5 - F3) < 0) {
  exact_xl <- -exact_xl
}

cat("\n")
cat("Ross index \n")
cat("beta distribution = ", formatC(beta_xl, digits = 4), "\n")
cat("   Exact solution = ", formatC(exact_xl, digits = 4), "\n")
    


leaf_area_data = read.csv('../Fortran/02_description_of_ecosystems/leaf_area/leaf_area_density.csv')

plot(leaf_area_data$LAD1, leaf_area_data$Height,  type = 'l', col = 'blue', xlab = 'Leaf area density (m^2 m^{-3})',
     ylab = 'Height (z/h_c)', main = 'Profiles')
lines(leaf_area_data$LAD2, leaf_area_data$Height, col = 'red')
lines(leaf_area_data$LAD3, leaf_area_data$Height, col = 'green')

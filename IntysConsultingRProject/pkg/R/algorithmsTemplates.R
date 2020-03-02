library(animation)
desc = c("This is a super cool example of Gradient Descent")
setwd("C:/Users/Gebruiker/Documents/Intys/RStudio/IntysConsultingRPRoject/OutputVisuals")
saveHTML({
  f1 = function(x, y) x^2 + 3 * sin(y)
  xx = grad.desc(f1, pi * c(-2, -2, 2, 2), c(-2 * pi, 2))
  
  xx$persp(col = "lightblue", theta = 30, phi = 30)
},
title = "Demo of Gradient Descent", description = desc, verbose = FALSE
)
save_gif({
  f1 = function(x, y) x^2 + 3 * sin(y)
  xx = grad.desc(f1, pi * c(-2, -2, 2, 2), c(-2 * pi, 2))
  
  xx$persp(col = "lightblue", theta = 30, phi = 30)
},
#,title = "Demo of Gradient Descent", description = desc, verbose = FALSE
)


#SOURCES
#https://www.rdocumentation.org/packages/gganimate/versions/1.0.5/topics/animate
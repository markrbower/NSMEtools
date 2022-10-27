library(ggplot2)
library(gganimate)
#> Loading required package: ggplot2

# We'll start with a static plot
p <- ggplot(iris, aes(x = Petal.Width, y = Petal.Length)) + 
  geom_point()

plot(p)

anim <- p + 
  transition_states(Species,
                    transition_length = 2,
                    state_length = 1)

anim


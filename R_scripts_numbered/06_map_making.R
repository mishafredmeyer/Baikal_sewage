library(ggmap)
library(grid)
library(viridis)


# 1. Load the data --------------------------------------------------------

metadata <- read.csv("../cleaned_data/metadata_20190320.csv", header = TRUE)

# 2. Create the map -------------------------------------------------------

map0 <- get_googlemap(c(108.171228, 53.538112), zoom = 6, maptype = "satellite")
map1 <- ggmap(map0) +
  geom_point(data = metadata, aes(long, lat), color = inferno(5)[3]) +
  #geom_path(data=xy, aes(x,y), color="red", lwd=1) +
  xlim(c(103, 110)) +
  ylim(c(51.5, 56)) +
  xlab("") +
  ylab("") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

map2 <- get_googlemap(c(105.199013, 51.944224), zoom = 10, maptype = "satellite")
map3 <- ggmap(map2) +
  geom_point(data = metadata, aes(x = long, y = lat), fill = inferno(5)[3],
             size = 10, stroke = 2, shape = 21) +
  xlim(c(104.8, 105.6)) +
  ylim(c(51.8, 52.1))


g1 <- ggplotGrob(map3)
grid.draw(g1)

pushViewport(viewport(x=0.25, y=0.8, w=.3, h=.3))
xy <- data.frame(x=c(105.4, 105.4, 105.6, 105.6, 105.4), 
                 y=c(51.8, 51.8, 51.9, 51.9, 51.8))
p2 <- ggmap(map1) + 
  geom_path(data=xy, aes(x,y), color="red", lwd=1) + 
  theme_void()
g2 <- ggplotGrob(map1)
grid.draw(g2)
grid.rect(gp=gpar(col="white", lwd=5))
popViewport()

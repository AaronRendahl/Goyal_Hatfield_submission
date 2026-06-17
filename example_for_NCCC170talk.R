library(tidyverse)
library(ggplot2)
library(patchwork)
set.seed(30)
rnormx <- function(n, sd) {
  sample(qnorm((1:n)/(n+1)))*sd
}
d <- expand_grid(rep=1:3, time=0:7) |> 
  mutate(m=4*exp(-.5*time), logm=log(m),
         e=rnormx(n(), 0.5), logy1=logm+e,
         logy=round(logy1*6)/6, y=exp(logy))
ggplot(d) + aes(time, m) + geom_point()
ggplot(d) + aes(time, logm) + geom_point()
ggplot(d) + aes(time, logy1) + geom_point()
ggplot(d) + aes(time, logy) + geom_point()

lm(logy ~ time, data=d)
# log(y) = 1.5 - 0.54 * time
# y = 4.5 exp (-0.54 * time)
p1 <- ggplot(d) + aes(time, y) + geom_point() + geom_smooth(method="lm", se=FALSE) +
  scale_y_log10(breaks=c(0.1, 0.2, 0.5, 1, 2, 5),
                minor_breaks=c(seq(0.01, 0.1, by=0.01),
                               seq(0.1, 1, by=0.1),
                               seq(1, 10, by=1)),
                labels=\(x) paste(x)) +
  labs(x="Time", y="Titer")

p2 <- ggplot(d) + aes(time, logy) + geom_point() + 
  geom_smooth(method="lm", se=FALSE) +
  scale_y_continuous(breaks=log(c(0.001, 1:5)), labels=0:5,
                     minor_breaks=log(seq(0.5, 10, by=0.5))) +
  expand_limits(y=log(0.001)) +
  coord_transform(y="exp") +
  labs(x="Time", y="Titer")

p2 | p1

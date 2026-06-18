library(tidyverse)
a <- readxl::read_excel("model_coefs.xlsx")
a |> summarize(minpwr=min(pwr), maxpwr=max(pwr), .by=Temp)
a |> summarize(minlrc=min(lrc), maxlrc=max(lrc), .by=Temp)
d |> distinct(Time, Temp)
#Temp4: 
t4x <- c(0:7, 0.5)
#TempRoom: 
trx <- c(0,2,4,6,8,10,12)

weibull <- function(x, Asym, Drop, lrc, pwr) {
  tibble(x=x, lrc=lrc, pwr=pwr, y=Asym-Drop*exp(-exp(lrc)*x^pwr))
}

x1 <- bind_rows(weibull(t4x, 0, -4, -4, 0.5),
          weibull(t4x, 0, -4, -4, 1),
          weibull(t4x, 0, -4, -4, 1.5))

x2 <- bind_rows(weibull(t4x, 0, -4, -2, 0.5),
                weibull(t4x, 0, -4, -2, 1),
                weibull(t4x, 0, -4, -2, 1.5))
xx <- bind_rows(x1, x2)

ggplot(xx) + aes(x, y, color=factor(pwr), shape=factor(lrc)) + 
  geom_point() +
  geom_line() +
  scale_y_log10()


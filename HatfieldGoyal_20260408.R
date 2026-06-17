# use preset starts for all Weibull fits and force > 0
# calculate AUC for each fit -- with a SE
# then do weighted ANOVA on AUC, on Subtype/Fomite, for each Temp separately

library(conflicted)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
conflicts_prefer(dplyr::filter)

d <- bind_rows(
  readxl::read_excel("AMPV_20260327.xlsx", 1) |> rename(Time="Days") |> mutate(Temp="4", .before=1),
  readxl::read_excel("AMPV_20260327.xlsx", 2) |> rename(Time="Hours") |> mutate(Temp="RT", .before=1)) |>
  pivot_longer(contains("_"), names_to=c("Subtype", "Rep"), names_sep="_", values_to='Titer6') |>
  mutate(Titer=Titer6/6) |>
  mutate(Group=case_when(Surface %in% c("AL", "PP", "PVC", "SS") ~ "Non-porous",
                         Surface %in% c("CB", "FB", "SF", "WO", "RB") ~ "Porous",
                         TRUE ~ "Standards") |>
           factor(levels=c("Standards", "Porous", "Non-porous")),
         .before="Surface") |>
  mutate(Surface=factor(Surface, levels=c("Media", "PBS", "AL", "PP", "PVC", "SS",
                                          "NL", "CB", "FB", "SF", "WO", "RB"))) |>
  filter(Time!=24)

weibull <- function(x, Asym, Drop, lrc, pwr) {
  Asym-Drop*exp(-exp(lrc)*x^pwr)  
}


# #####
# ## how to do just one
# 
# do_temp <- "RT"
# idx <- 23
# d1 <- d |> filter(Temp==do_temp)
# dxs <- d1 |> nest(.by=c(Temp, Surface, Subtype, Group))
# 
# 
# dx <- dxs$data[[idx]]
# 
# mx <- nls(Titer ~ weibull(Time, Asym, Drop, lrc, pwr), data=dx,
#           start = c(Asym=0, Drop=-2, lrc=-2, pwr=1),
#           lower=c(0,-Inf, -Inf, -Inf), upper=c(Inf, 0, Inf, Inf),
#           algorithm="port")
# 
# px <- profile(mx)
# par(mfrow=c(2,2)); plot(px)
# 
# fx <- tibble(Time=sort(unique(dx$Time)))
# fx$Titer <- predict(mx, newdata=fx)
# ggplot(dx) + aes(Time, Titer) + geom_point() +
#   geom_line(data=fx) 
# 
# broom::tidy(mx) |> select(term, est=estimate, se=std.error)
# vcov(mx)

#### just get fits #################
fitWeibull <- function(dx, Temp, start) {
  if(!missing(Temp)) {
    if(Temp=="4") {
      start <- c(Asym=0, Drop=-4, lrc=-2, pwr=1)
    } else {
      start <- c(Asym=0, Drop=-(max(dx$Titer)-min(dx$Titer)), lrc=-2, pwr=1)
    }
  }
  nls(Titer ~ weibull(Time, Asym, Drop, lrc, pwr), data=dx,
            start = start,
            lower=c(0,-Inf, -Inf, -Inf), upper=c(Inf, 0, Inf, Inf),
            algorithm="port")
}

d |> nest(.by=c(Temp, Surface, Subtype, Group)) |> 
  mutate(model=map2(data, Temp, fitWeibull),
         coef=map(model, \(x) as_tibble_row(coef(x)))) |>
  unnest(coef) |>
  select(-c(data, model)) |>
  writexl::write_xlsx("model_coefs.xlsx")

#### do simulation ##########################
simAUC <- function(dx, show=FALSE, N=2000,
                   #start=c(Asym=0, Drop=-4, lrc=-2, pwr=1) # for 4
                   start=c(Asym=0, Drop=-(max(dx$Titer)-min(dx$Titer)), lrc=-2, pwr=1) # for RT
                   ) {
  mx <- nls(Titer ~ weibull(Time, Asym, Drop, lrc, pwr), data=dx,
            start = start,
            lower=c(0,-Inf, -Inf, -Inf), upper=c(Inf, 0, Inf, Inf),
            algorithm="port")
  
  fx <- tibble(Time=sort(unique(dx$Time)))
  fx$Titer <- predict(mx, newdata=fx)
  ax <- fx |> arrange(Time) |>
    mutate(my=(Titer+lead(Titer))/2, dx=lead(Time)-Time, a = my*dx) |>
    summarize(AUC=sum(a, na.rm=TRUE)) |>
    bind_cols(fx |> filter(Time==0) |> select(Titer0=Titer)) |>
    mutate(aAUC=AUC/Titer0)
  
  dx$fit <- fitted(mx)
  sigma <- summary(mx)$sigma
  
  to_titer <- \(x) pmax((round((x*6+1)/2)*2-1)/6, 0L)  
  cx <- tibble(term=names(coef(mx)))
  for(k in 1:N) {
    #if(k %% 100 == 0) message(k)
    mxx <- 0L
    ## trouble here is that every so often, get a simulation that can't run
    ## let's just skip those
    while(isTRUE(all.equal(mxx, 0L))) {
      mxx <- tryCatch({
        dx$sim <- rnorm(nrow(dx), dx$fit, sigma) |> to_titer()
        mxx <- nls(sim ~ weibull(Time, Asym, Drop, lrc, pwr), data=dx,
                   start = c(Asym=0, Drop=-4, lrc=-2, pwr=1),
                   lower=c(0,-Inf, -Inf, -Inf), upper=c(Inf, 0, Inf, Inf),
                   algorithm="port")
      }, error=\(e) 0L)
    }
    fx[[sprintf("sim_%03d",k)]] <- predict(mxx, newdata=fx)
    cx[[sprintf("sim_%03d",k)]] <- coef(mxx)
    
    # Drop0 = -mean(dx$sim[dx$Time==0])
    # mxx <- nls(sim ~ weibull(Time, Asym, Drop, lrc, pwr), data=dx,
    #            start = c(Asym=0, Drop=Drop0, lrc=-3, pwr=2.5),
    #            lower=c(0,-Inf, -Inf, -Inf), upper=c(Inf, 0, Inf, Inf),
    #            algorithm="port")
    # 
    # fx <- tibble(Time=sort(unique(dx$Time)))
    # fx$sim <- weibull(fx$Time, 0, Drop0, -3, 2.5)
    # ggplot(dx) + aes(Time, sim) + geom_point() +
    #   geom_line(data=fx)
  }
  message("done")
  
  fxs <- fx |> select(Time, starts_with("sim")) |>
    pivot_longer(starts_with("sim"), values_to="Titer", names_to="Sim")
  

  px <- ggplot(dx) + aes(Time, Titer) +
    geom_line(data=fxs, aes(group=Sim), alpha=0.1) +
    geom_point() +
    geom_line(data=fx)
  if(show) {
    plot(px)
  }
  axs <- fxs |> arrange(Sim, Time) |>
    mutate(my=(Titer+lead(Titer))/2, dx=lead(Time)-Time, a = my*dx, .by=Sim) |>
    summarize(AUC=sum(a, na.rm=TRUE), .by=Sim) |> 
    left_join(fxs |> filter(Time==0) |> select(Sim, Titer0=Titer), by="Sim") |>
    mutate(aAUC=AUC/Titer0)
  tibble_row(model=list(mx), plot=list(px), fit=list(fx), AUC=ax$aAUC, mAUC=mean(axs$aAUC), seAUC=sd(axs$aAUC))
}

for(do_temp in c("RT", "4")) {
  
temp_code <- if(do_temp=="4") "T4" else "TR"
temp_text <- if(do_temp=="4") "Temp 4" else "Room Temp"

d1 <- d |> filter(Temp==do_temp)
dxs <- d1 |> nest(.by=c(Temp, Surface, Subtype, Group))

simfile <- sprintf("sim_%s.RDS", temp_code)
if(!file.exists(simfile)) {
  fits <- dxs |> mutate(map_dfr(data, simAUC))
  write_rds(fits, file=sprintf("sim_%s.RDS", temp_code))
}
fits <- read_rds(simfile)

f1s <- fits |> select(Temp, Surface, Subtype, Group, fit) |> unnest(fit) |> select(-Titer) |>
  pivot_longer(starts_with("sim_"), names_to='Sim', values_to='Titer') |>
  mutate(Sim=as.integer(str_remove(Sim, "sim_"))) |>
  mutate(Titer0=Titer[Time==0], TiterNorm=Titer/Titer0,
         .by=c(Temp, Surface, Subtype, Group, Sim)) |>
  mutate(G=sprintf("%s_%s_%s_%03d", Temp, Surface, Subtype, Sim) |> factor(), .after='Sim')
f1 <- fits |> select(Temp, Surface, Subtype, Group, fit) |> unnest(fit) |> select(-starts_with("sim")) |>
  mutate(Titer0=Titer[Time==0], TiterNorm=Titer/Titer0,
         .by=c(Temp, Surface, Subtype, Group))
f1ci <- f1s |> summarize(Titer_lwr=quantile(Titer, 0.025), Titer_upr=quantile(Titer, 0.975),
                         TiterNorm_lwr=quantile(TiterNorm, 0.025), TiterNorm_upr=quantile(TiterNorm, 0.975),
                         .by=c(Temp, Surface, Subtype, Group, Time))

# ggplot(d1) + aes(Time, Titer) + facet_grid(Subtype ~ Surface) +
#   geom_line(data=f1s, aes(group=G), alpha=0.05) +
#   geom_point() +
#   geom_line(data=f1)

xb <- sort(unique(d1$Time))[c(TRUE, FALSE)]

px1 <- ggplot() + aes(Time, Titer, color=Subtype) + 
  facet_wrap( ~ Surface, nrow=2) + 
  #geom_line(data=f1s, aes(group=G), alpha=0.1) +
  geom_ribbon(data=f1ci, aes(y=NULL, ymin=Titer_lwr, ymax=Titer_upr, fill=Subtype), linewidth=0, alpha=0.5) +
  geom_line(data=f1) +
  geom_point(data=d1, fill='white', pch=21, group=1, 
             position=ggbeeswarm::position_beeswarm()) +
  scale_x_continuous(breaks=xb) +
  labs(title=sprintf("Fitted Weibull curves, %s", temp_text)) +
  theme(legend.position = "bottom") +
  labs(y="Titer")
ggsave(px1, filename=sprintf("%s_plot1_Fits.pdf", temp_code), width=12, height=5)

px2 <- ggplot() + aes(Time, TiterNorm, color=Subtype) + 
  facet_wrap( ~ Surface, nrow=2) +
  #geom_line(data=f1s, aes(group=G), alpha=0.1) +
  geom_ribbon(data=f1ci, aes(y=NULL, ymin=TiterNorm_lwr, ymax=TiterNorm_upr, fill=Subtype), linewidth=0, alpha=0.5) +
  geom_line(data=f1) +
  annotate(x=0.2, y=0.4, label="AUC (\u00B1 SE)", hjust=0, size=10/.pt, geom="text") +
  geom_text(aes(label=sprintf("%0.2f \u00B1 %0.2f", AUC, seAUC), x=0.2, y=0.2+(Subtype=='A')*0.1), data=fits,
            hjust=0, size=10/.pt) +
  scale_y_continuous(breaks=seq(0,1,by=0.2), labels=scales::label_percent()) +
  scale_x_continuous(breaks=xb) +
  labs(title=sprintf("Rescaled Weibull curves, %s", temp_text),
       y="Normalized Titer (to Time 0)") +
  theme(legend.position = "bottom")
ggsave(px2, filename=sprintf("%s_plot2_NormAUC.pdf", temp_code), width=12, height=5)

z95 <- qnorm(1-0.05/2)
px3 <- ggplot(fits) + aes(Surface, AUC, color=Subtype) + 
  facet_grid(.~Group, scales="free", space="free") +
  geom_pointrange(aes(y=AUC, ymin=AUC-z95*seAUC, ymax=AUC+z95*seAUC),
                  position=position_dodge(width=0.3)) +
  ylab("AUC (95% CI)") +
  labs(title=temp_text)  +
  theme(legend.position = "bottom")
ggsave(px3, filename=sprintf("%s_plot3_AUC.pdf", temp_code), width=6, height=4)


ztest <- function(m1, m2, s1, s2) {
  s <- sqrt(s1^2+s2^2)
  d <- m2-m1
  z <- d/s
  p <- pnorm(abs(z), lower.tail=FALSE)
  tibble(diff=d, diff_se=s, z=z, p.value=p)
}

getpairwise <- function(ff, n="results") {
  zp2 <- t(combn(as.character(ff$ID), 2)) |> 
    as.data.frame() |> as_tibble() |>
    mutate(map2_dfr(V1, V2, \(v1, v2) {
      f1 <- ff |> filter(ID==v1)
      f2 <- ff |> filter(ID==v2)
      ztest(f1$AUC, f2$AUC, f1$seAUC, f2$seAUC)
    })) |>
    mutate(p.adjust=p.adjust(p.value))
  
  xx <- multcompView::multcompLetters(setNames(zp2$p.adjust, paste(zp2$V1, zp2$V2, sep="-")))
  zp1 <- ff |> left_join(tibble(ID=names(xx$monospacedLetters), .group=unname(xx$monospacedLetters)), by="ID")
  out <- list(AUC=zp1, pairwise=zp2)
  names(out) <- paste(n,c("AUC", "pairwise"), sep="-")
  out
}

day3 <- c("NL", "PBS")

ffA <- fits |> select(Temp, Surface, Subtype, Group, AUC, seAUC) |>
  arrange(AUC) |>
  mutate(ID=paste(Temp, Surface, Subtype, sep=".")) |>
  filter(Subtype=="A", !Surface %in% day3) |>
  getpairwise("A")

ffB <- fits |> select(Temp, Surface, Subtype, Group, AUC, seAUC) |>
  arrange(AUC) |>
  mutate(ID=paste(Temp, Surface, Subtype, sep=".")) |>
  filter(Subtype=="B", !Surface %in% day3) |>
  getpairwise("B")

ff3 <- fits |> select(Temp, Surface, Subtype, Group, AUC, seAUC) |>
  arrange(AUC) |>
  mutate(ID=paste(Temp, Surface, Subtype, sep="."), .before=1) |>
  filter(Surface %in% day3) |>
  getpairwise("Day3")

writexl::write_xlsx(c(ffA, ffB, ff3),
                    path=sprintf("%s_AUC.xlsx", temp_code))

fitsx <- fits |> filter(!Surface %in% day3)
fx <- bind_rows(ffA$`A-AUC`, ffB$`B-AUC`)
px4 <- ggplot(fitsx) + aes(Surface, AUC, color=Group) + 
  facet_wrap(.~Subtype) +
  geom_pointrange(aes(y=AUC, ymin=AUC-z95*seAUC, ymax=AUC+z95*seAUC),
                  position=position_dodge(width=0.3)) +
  ylab("Normalized AUC (95% CI)") +
  labs(title=temp_text)  +
  theme(legend.position = "bottom") +
  scale_x_discrete(limits=rev) +
  scale_y_continuous(breaks=xb) +
  coord_flip()
px4txt <- px4 + geom_text(aes(label=.group), y=Inf, data=fx, color="black", size=8/.pt, hjust=1)
ggsave(px4, filename=sprintf("%s_plot4x_AUC_compare.pdf", temp_code), width=6, height=4)
ggsave(px4txt, filename=sprintf("%s_plot4_AUC_compare.pdf", temp_code), width=6, height=4)

fitsx <- fits |> filter(Surface %in% day3) |>
  mutate(X=paste(Surface, Subtype))
fx <- ff3$`Day3-AUC` |>
  mutate(X=paste(Surface, Subtype))
px4b <- ggplot(fitsx) + aes(X, AUC, color=Subtype) + 
  geom_pointrange(aes(y=AUC, ymin=AUC-z95*seAUC, ymax=AUC+z95*seAUC),
                  position=position_dodge(width=0.3)) +
  ylab("Normalized AUC (95% CI)") +
  labs(title=temp_text)  +
  theme(legend.position = "bottom") +
  scale_x_discrete(limits=rev, name="Surface / Subtype") +
  coord_flip()
px4btxt <-px4b + geom_text(aes(label=.group), y=Inf, data=fx, color="black", size=8/.pt, hjust=1)
ggsave(px4b, filename=sprintf("%s_plot4bx_AUC_compare.pdf", temp_code), width=4, height=3)
ggsave(px4btxt, filename=sprintf("%s_plot4b_AUC_compare.pdf", temp_code), width=4, height=3)

px5 <- ggplot(f1) + aes(Time, TiterNorm, color=Surface) +
  facet_grid(. ~ Subtype) +
  geom_line() +
  labs(title=sprintf("Rescaled Weibull curves, %s", temp_text)) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks=xb)

ggsave(px5, filename=sprintf("%s_plot5_AUC_overlap.pdf", temp_code), width=6, height=4)

}
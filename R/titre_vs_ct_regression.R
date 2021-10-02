library(data.table)
library(readxl)
library(tidyverse)

LOD = 20

#### Read data in and clean ####

file='J, K, L, M Master Logbook_V4.xlsx'

jklm <- list(
  read_excel(file, sheet = 'J Set'),
  read_excel(file, sheet = 'K Set'),
  read_excel(file, sheet = 'L Set'),
  read_excel(file, sheet = 'M Set')
)

joined_char <- lapply(seq_len(length(jklm)), function(i) {
  jklm[[i]] %>%
    select(-Date) %>%
    mutate_all(as.character) %>%
    ## M set were in a different batch
    mutate(set=case_when(i<4 ~ 0, TRUE ~ 1))
}) %>%
  bind_rows() %>%
  mutate_all(function(v) replace(v, which(trimws(toupper(v))=='NDET'), 45)) %>%
  mutate_all(function(v) replace(v, which(trimws(toupper(v))=='NA'), NA))  %>%
  mutate_all(function(v) replace(v, which(trimws(toupper(v))=='UNDETERMINED'), NA)) 

## check fixed all the non-numeric values
all(is.na(lapply(joined_char[,c('E Ct', 'sgE Ct', 'Direct Titer (Veros)', 
                      'Direct Titer (TMPRSS2)', 'Direct Titer (Vero)')], 
       function(col_vals) {
         unique(col_vals[which(is.na(as.numeric(col_vals)))])
       })))

joined <- cbind(
  joined_char[,c('CID', 'Code', 'Variant', 'set')],
  joined_char[,c('E Ct', 'sgE Ct', 'Direct Titer (TMPRSS2)')] %>%
    mutate_all(as.numeric)
) %>% as.data.table() %>%
  (function(df) {
  df[
    `E Ct` <30, #remove lone high-CT point that wasn't planned to be titred
  ][gsub('\\.', '', Variant) %in% c('B117', 'B1429', 'B16172', 'B161722') & 
           !is.na(`Direct Titer (TMPRSS2)`)
  ][
    ,Variant := factor(
      Variant, 
      levels = c('B.1.1.7', 'B.1.429', 'B.1.617.2', 'B.1.617.2.2'),
      labels = c('Alpha', 'Epsilon', 'Delta', 'Delta')
    )
  ]})()


#### Modeling ####
####### gaussian ####
#### CT E functions ####
fit_model_e <- function(df) {
  model <- df[,{
    model=glm(log10(`Direct Titer (TMPRSS2)`) ~ `E Ct` + Variant, data = .SD, family='gaussian')
  }]
  print(summary(model))
  cat('Original scale coefficients:\n')
  print(t(t(10^model$coefficients)))
  return(model)
}

predict_and_join <- function(df, model) {
  ## function to predict and tack on those predictions to a data frame where model was possibly
  ## fitted using different data
  df[,{
    fitted = predict.glm(model, .SD, se.fit=TRUE)
    link_inverse = function(x) 10^x
    c(.SD, list(predicted = link_inverse(fitted$fit), 
                
                LB=link_inverse(fitted$fit - qnorm(0.975) * fitted$se.fit), 
                UB =link_inverse(fitted$fit + qnorm(0.975) * fitted$se.fit)))
  }]
}

### trying various different models to see how variable coefficient estimates are
### full_lod estimate
full_lod_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD]
full_lod_model <- fit_model_e(full_lod_data)

#sensitivity to missing value selection
invisible(fit_model_e(copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=0.1]))
invisible(fit_model_e(copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=1]))
invisible(fit_model_e(copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=10]))

invisible(fit_model_e(copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD]))
invisible(fit_model_e(copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD/10]))
invisible(fit_model_e(copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD/100]))


half_lod_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD/2]
half_lod_model <- fit_model_e(half_lod_data)
glm(log10(`Direct Titer (TMPRSS2)`) ~ `E Ct` + Variant + set, data = half_lod_data, family='gaussian') %>% summary()
glm(log10(`Direct Titer (TMPRSS2)`) ~ `E Ct` * Variant, data = half_lod_data, family='gaussian') %>% summary()
glm(`Direct Titer (TMPRSS2)` ~ `E Ct` + Variant, data = joined, family='quasipoisson') %>% summary()

### low (LOD/10) is the reference model
low_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD/10]
low_model <- fit_model_e(low_data)



set.seed(42)
randomized_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=10^runif(.N, log10(1.5), log10(2.5))]
# randomized_model <- fit_model_e(randomized_data)

trimmed_data <- copy(joined)[`Direct Titer (TMPRSS2)`!=0,]
trimmed_model <- fit_model_e(trimmed_data)


ggplot(predict_and_join(randomized_data, low_model), 
       aes(x=`E Ct`, y=`Direct Titer (TMPRSS2)`, color=Variant,
           shape = Variant,
           group = Variant)) +
  geom_ribbon(aes(ymin=LB, ymax=UB, fill=Variant, color=NULL), alpha=0.14 ) +
  geom_line(aes(y=predicted)) +
  geom_point(size=2.5, alpha=0.9) +
  geom_hline(yintercept = LOD, linetype=2) +
  scale_size_identity() +
  scale_color_manual(values=rcartocolor::carto_pal(name='Prism')[c(3,7,1)]) +
  scale_fill_manual(values=rcartocolor::carto_pal(name='Prism')[c(3,7,1)]) +
  theme_bw() +
  scale_y_log10(breaks=10^c(0:6)) +
  coord_cartesian(ylim=c(1, 1e6))
ggsave(paste0('lm_e_', Sys.Date(), '.pdf'), height = 8, width = 10)



#### CT sgE ####
# functions:
fit_model_sge <- function(df) {
  model <- df[,{
    model=glm(log10(`Direct Titer (TMPRSS2)`) ~ `sgE Ct` + Variant, data = .SD, family='gaussian')
  }]
  print(summary(model))
  cat('Original scale coefficients:\n')
  print(t(t(10^model$coefficients)))
  return(model)
}

### trying various different models to see how variable coefficient estimates are
### full_lod estimate
full_lod_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD]
full_lod_model <- fit_model_sge(full_lod_data)

half_lod_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD/2]
half_lod_model <- fit_model_sge(half_lod_data)
glm(log10(`Direct Titer (TMPRSS2)`) ~ `sgE Ct` * Variant, data = half_lod_data, family='gaussian') %>% summary()
glm(log10(`Direct Titer (TMPRSS2)`) ~ `sgE Ct` + Variant +set, data = half_lod_data, family='gaussian') %>% summary()

## ref model
low_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=LOD/10]
low_model <- fit_model_sge(low_data)

set.seed(42)
randomized_data <- copy(joined)[`Direct Titer (TMPRSS2)`==0,`Direct Titer (TMPRSS2)`:=10^runif(.N, log10(1.5), log10(2.5))]
randomized_model <- fit_model_sge(randomized_data)

trimmed_data <- copy(joined)[`Direct Titer (TMPRSS2)`!=0,]
trimmed_model <- fit_model_sge(trimmed_data)


### plot ###
ggplot(predict_and_join(randomized_data, low_model), 
       aes(x=`sgE Ct`, y=`Direct Titer (TMPRSS2)`, color=Variant, 
           shape = Variant, 
           group = Variant)) +
  geom_ribbon(aes(ymin=LB, ymax=UB, fill=Variant, color=NULL), alpha=0.14 ) +
  geom_line(aes(y=predicted)) +
  geom_point(size=2.5, alpha=0.9) +
  geom_hline(yintercept = LOD, linetype=2) +
  scale_size_identity() +
  scale_color_manual(values=rcartocolor::carto_pal(name='Prism')[c(3,7,1)]) +
  scale_fill_manual(values=rcartocolor::carto_pal(name='Prism')[c(3,7,1)]) +
  theme_bw() +
  scale_y_log10(breaks=10^c(0:6)) +
  coord_cartesian(ylim=c(1, 1e6))
ggsave(paste0('lm_sgE_', Sys.Date(), '.pdf'), height = 8, width = 10)


#Illustrate Z measure
options.posterior <- seq(0,1,0.01)
results05 <- Z(options.posterior,0.5)
results02 <- Z(options.posterior,0.2)
results08 <- Z(options.posterior,0.8)
p05<- ggplot()+
  geom_line(aes(x=options.posterior,y=results05))+theme_minimal(base_size = 6)+
  labs(title="Z as a function of posterior",subtitle="With prior = 0.5")
p02<- ggplot()+
  geom_line(aes(x=options.posterior,y=results02))+theme_minimal(base_size = 6)+
  labs(title="Z as a function of posterior",subtitle="With prior = 0.2")
p08<- ggplot()+
  geom_line(aes(x=options.posterior,y=results08))+theme_minimal(base_size = 6)+
  labs(title="Z as a function of posterior",subtitle="With prior = 0.8")
grid.arrange(p02,p05,p08,ncol=3)

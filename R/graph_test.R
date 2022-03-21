library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(magrittr)
pdf(file = "plots/lcpa.pdf", width = 4, height = 6)

grViz("
digraph {
# latent variables
lprofile [shape = oval]

#regressions
lprofile -> lclass1
lprofile -> lclass2
lprofile -> lclass3

lclass1 [shape = oval];
lclass2 [shape = oval];
lclass3 [shape = oval];

item11 [shape=box];
item21 [shape=box];
item31 [shape=box];

item12 [shape=box];
item22 [shape=box];
item32 [shape=box];

item13 [shape=box];
item23 [shape=box];
item33 [shape=box];

lclass1 -> item11;
lclass1 -> item21;
lclass1 -> item31;

lclass2 -> item12;
lclass2 -> item22;
lclass2 -> item32;

lclass3 -> item13;
lclass3 -> item23;
lclass3 -> item33;

}
") %>% export_svg %>% charToRaw %>% rsvg_pdf("plots/lcpa2.pdf")


grViz("
digraph {
# latent variables
lcluster [shape = oval]

#regressions
lcluster -> lclass

#measurement model for latent a
subgraph cluster0{
label = 'in each group'
labeljust='r';

lclass [shape = oval];
item1 [shape=box];
item2 [shape=box];
item3 [shape=box];
}


lclass -> item1;
lclass -> item2;
lclass -> item3;

}
") %>% export_svg %>% charToRaw %>% rsvg_pdf("plots/mlca.pdf")

grViz("
digraph {
# latent variables
lcluster0 [shape = oval]

#regressions
lcluster0 -> lcluster1
lcluster1 -> lclass

#measurement model for latent a
subgraph cluster0{
   label = 'in each school';
   labeljust='r';
   lcluster1 [shape = oval]
   subgraph cluster1{
      label = 'in each class';
      labeljust='r';

      lclass [shape = oval];
      item1 [shape=box];
      item2 [shape=box];
      item3 [shape=box];
   }
}

lclass -> item1;
lclass -> item2;
lclass -> item3;

}
") %>% export_svg %>% charToRaw %>% rsvg_pdf("plots/mlca3.pdf")



grViz("
digraph {
# latent variables
lcluster [shape = oval]

#regressions
lcluster -> lprofile

#measurement model for latent a
subgraph cluster0{
label = 'in each group'
labeljust = 'r';
# latent variables
lprofile [shape = oval]

#regressions
lprofile -> lclass1
lprofile -> lclass2
lprofile -> lclass3

lclass1 [shape = oval];
lclass2 [shape = oval];
lclass3 [shape = oval];

item11 [shape=box];
item21 [shape=box];
item31 [shape=box];

item12 [shape=box];
item22 [shape=box];
item32 [shape=box];

item13 [shape=box];
item23 [shape=box];
item33 [shape=box];
}



lclass1 -> item11;
lclass1 -> item21;
lclass1 -> item31;

lclass2 -> item12;
lclass2 -> item22;
lclass2 -> item32;

lclass3 -> item13;
lclass3 -> item23;
lclass3 -> item33;

}
") %>% export_svg %>% charToRaw %>% rsvg_pdf("plots/mlcpa.pdf")


grViz("
digraph {
# latent variables
lgroup [shape = oval]

#regressions
lgroup -> lclass

# #measurement model for latent a
# subgraph cluster0{
# label = D
# labeljust='r';

lclass [shape = oval];

item1 [shape=box];
item2 [shape=box];
item3 [shape=box];

gitem1 [shape=box];
gitem2 [shape=box];
gitem3 [shape=box];

# }


lclass -> item1;
lclass -> item2;
lclass -> item3;

gitem1 -> lgroup [dir = back];
gitem2 -> lgroup [dir = back];
gitem3 -> lgroup [dir = back];
}
") %>% export_svg %>% charToRaw %>% rsvg_pdf("plots/lcalg.pdf")



grViz("
digraph {
# latent variables
lgroup[shape = oval]
lprofile [shape = oval]

#regressions
lgroup -> lprofile
lprofile -> lclass1
lprofile -> lclass2
lprofile -> lclass3

# #measurement model for latent a
# subgraph cluster0{
# label = D
# labeljust='r';

lclass1 [shape = oval];
lclass2 [shape = oval];
lclass3 [shape = oval];

item11 [shape=box];
item21 [shape=box];
item31 [shape=box];

item12 [shape=box];
item22 [shape=box];
item32 [shape=box];

item13 [shape=box];
item23 [shape=box];
item33 [shape=box];

item1 [shape=box];
item2 [shape=box];
item3 [shape=box];
# }


lclass1 -> item11;
lclass1 -> item21;
lclass1 -> item31;

lclass2 -> item12;
lclass2 -> item22;
lclass2 -> item32;

lclass3 -> item13;
lclass3 -> item23;
lclass3 -> item33;

item1 -> lgroup [dir = back];
item2 -> lgroup [dir = back];
item3 -> lgroup [dir = back];
}
") %>% export_svg %>% charToRaw %>% rsvg_pdf("plots/lcpalg.pdf")

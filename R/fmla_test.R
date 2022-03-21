fy = list(
   l1(3) ~ y11 + y21 + y31 + y41,
   l2(3) ~ y12 + y22 + y32 + y42,
   l3(3) ~ y13 + y23 + y33 + y43,
   pf(3) ~ l1 + l2 + l3,
)

fx = list(pf ~ x1 + x2)

fx = list(pf ~ x + (1 + z | grp)
)

const = list(
   eq(y11, y12, y13),
   eq(y21, y22, y23),
   eq(y31, y32, y33)
)

catlvm(fy, fx, const, data = dataset)


####################################

fy = list(
   lclass("lclss", 3, grp, "inv") ~ y11 + y21 + y31 + y41,
   lclass("lclst", 3) ~ lclss
)

fx = list(
   lclass ~ x + (z | grp)
)

const = list(
   y11 ~ y12 ~ y13,
   y21 ~ y22 ~ y23,
   y31 ~ y32 ~ y33
)

catlvm(fy, fx, const, data = dataset)


##########################################
##########################################

lv = define_lv(
   c("l1", 3), c("l2", 3), c("l3", 3), c("p", 3)
)

path = define_path(
   l1 %by% c(y11, y21, y31, y41),
   l2 %by% c(y12, y22, y32, y42),
   l3 %by% c(y13, y23, y33, y43),
   p %by% c(l1, l2, l3),
   p ~ x1 + x2 + x3
)

paths = list(
   l2o(from = "lc1", to = c("y11", "y21", "y31", "y41")),
   l2o(from = "lc2", to = c("y12", "y22", "y32", "y42")),
   l2o(from = "lc3", to = c("y13", "y23", "y33", "y43")),
   l2l(from = "pf", to = c("lc1", "lc2", "lc3")),
   x2l(from = c("x1", "x2", "x3"), to = c("pf"))
)

catlvm(lv, path, data = dataset)

##########################################

observed = list(
   lc1(3) ~ y11 + y21 + y31 + y41,
   lc2(3) ~ y12 + y22 + y32 + y42,
   lc3(3) ~ y13 + y23 + y33 + y43
)

latent = list(
   pf(3) ~ lc1 + lc2 + lc3
)

covariate = list(
   pf ~ x + (1 + z | grp)
)

catlvm(observed, latent, covariate, data = dataset)


# The R and Python packages give the same numbers on the bundled mock city.
# Reference values come from the Python package rd2dnet 0.1.0.

test_that("coordinates match the Python package", {
  expect_equal(mock_coords()$s[1:6], c(829.087, 1649.546, -1353.994, -1783.17, 1648.397, 1653.963),
               tolerance = 1e-12)
})

test_that("estimates and standard errors match the Python package", {
  y <- mock_city()$obs$outcome
  cc <- mock_coords()
  res <- rdnet(y, cc, h = 600, vce = "hc1", seed = 0)
  ok <- which(res$table$status == "ok")[1:4]
  expect_equal(ok, 3:6)
  t <- res$table[ok, ]
  expect_equal(t$tau, c(2.2879700352572083, 2.226131870704858, 2.2374576858442436, 2.216318573301832),
               tolerance = 1e-10)
  expect_equal(t$se, c(0.6423910522468664, 0.519468183514787, 0.3817868050575414, 0.450498262570776),
               tolerance = 1e-10)
  expect_equal(t$tau_bc, c(2.5630875033719818, 2.671419704490127, 2.338716529414876, 2.3138693971233018),
               tolerance = 1e-10)
  expect_equal(t$se_rb, c(0.9542401673002789, 0.6640019263539401, 0.6604308274435798, 0.7795185876927633),
               tolerance = 1e-10)
  conley <- rdnet(y, cc, h = 600, conley_cutoff = 150, seed = 0)$table$se_rb[ok]
  expect_equal(conley, c(0.8979338736065969, 0.6601068377141024, 0.5275681466518771, 0.6208458892927217),
               tolerance = 1e-10)
  cluster <- rdnet(y, cc, h = 600, cluster = mock_city()$obs$block, seed = 0)$table$se_rb[ok]
  expect_equal(cluster, c(0.3221750701567485, 0.4308191219209627, 0.45480015406506874, 0.45409127249969544),
               tolerance = 1e-10)
})

test_that("the density test matches the Python package", {
  r <- density(mock_coords(), normalize = "network")
  expect_equal(c(r$diff_bc, r$se, r$h_left), c(3.370389666344076e-06, 3.200746336495476e-06, 12249.64743015767),
               tolerance = 1e-8)
  r <- density(mock_coords(), normalize = "none")
  expect_equal(c(r$diff_bc, r$se, r$h_left), c(-5.152889755147277e-05, 9.395949397424248e-05, 326.7653651609291),
               tolerance = 1e-8)
})

hhg.example.datagen = function(n, example) {
  if (example == '') {
  } else if (example == '4indclouds') {
    .datagen4indclouds(n)
  } else if (example == '2Parabolas') {
    .datagen2Parabolas(n)
  } else if (example == 'W') {
    .datagenW(n)
  } else if (example == 'Parabola') {
    .datagenParabola(n)
  } else if (example == 'Diamond') {
    .datagenDiamond(n)
  } else if (example == 'Circle') {
    .datagenCircle(n)
  } else {
    stop('Unexpected example specified. Please consult the documentation.')
  }
}

.datagen4indclouds = function(n) {
  dx = rnorm(n) / 3
  dy = rnorm(n) / 3
  cx = sample(c(-1, 1), size = n, replace = T)
  cy = sample(c(-1, 1), size = n, replace = T)
  u = cx + dx
  v = cy + dy 
  return (rbind(u, v))
}

.datagen2Parabolas = function(n) {
  x = seq(-1, 1, length = n)
  y = (x ^ 2 + runif(n) / 2) * (sample(c(-1, 1), size = n, replace = T))
  return (rbind(x, y))
}

.datagenW = function(n) {
  x = seq(-1, 1, length = n)
  u = x + runif(n)/3
  v =  4*( ( x^2 - 1/2 )^2 + runif(n)/500 )
  return (rbind(u,v))
}

.datagenParabola = function(n) {
  x = seq(-1, 1, length = n)
  y = (x ^ 2 + runif(n)) / 2
  return (rbind(x,y))
}

.datagenDiamond = function(n) {
  x = runif(n, min = -1, max = 1)
  y = runif(n, min = -1, max = 1)

  theta = -pi / 4
  rr = rbind(c(cos(theta), -sin(theta)),
             c(sin(theta),  cos(theta)))
  tmp = cbind(x, y) %*% rr
  u = tmp[,1]
  v =  tmp[,2]
  return (rbind(u, v))
}

.datagenCircle = function(n) {
  x = seq(-1, 1, length = n)
  u = sin(x * pi) + rnorm(n) / 8
  v = cos(x * pi) + rnorm(n) / 8
  return (rbind(u, v))
}

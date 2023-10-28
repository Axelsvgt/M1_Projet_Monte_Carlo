#################################
###### PROJET MONTE CARLO #######
#################################

# Adrien PASSUELLO & Axel SAUVAGET
# Groupe : 117


# NB : delta doit être environ égal à 7%
set.seed(243)# on définit une graine 
# (attention sur le document PDF les résultats peuvent différer
# alors que la graine est la même)

####### EXERCICE 1 #######


###### Question 1 ######
# Densité f d'une loi gamma de paramètre de forme alpha = 3 et d'échelle theta = 2
densite_f <- function(x,alpha = 3,theta = 2){
  return(dgamma(x, shape = alpha, scale = theta))
}

# Densité g de paramètre mu = -13 et b = 2
densite_g <- function(y, mu = -13, b = 2){
  return((y-mu+2*b)/(4*b*b) * (y >= mu-2*b) * (y <= mu) + exp(-(y-mu)/b)/(2*b) * (y > mu))
}

# Fonction de répartition de G calculée théoriquement
fdr_G <- function(y,mu = -13, b = 2){
  return((y - (mu-2*b))^2/(8*b*b)*(mu - 2* b <= mu)*(y <= mu)  +  (1-exp(-(y-mu)/b)/2)*(y>mu))
}

# Fonction inverse généralisée de G calculée théoriquement
G_fleche <- function(u,mu = -13, b = 2){
  return((2*b*(sqrt(2*u) - 1) + mu)*(0 <= u)*(u <= 1/2)  +  (mu-b*log(2*(1-u)))*(u > 1/2)*(u < 1) )
}

# Fonction de simulation de n v.a. de densité g (via la la méthode de la fonction inverse)
rgen_g <- function(n, mu = -13, b = 2){
  # On simule n v.a. uniforme sur [0,1]
  U <- runif(n)
  # On retourne G_flèche(U) car si U ~ U([0,1]), alors G_flèche(U) et Y sont de même loi
  return(G_fleche(U,mu,b))
}

# Fonction de simulation de n v.a. de densité f (via la fonction de R rgamma)
rgen_f <- function(n, alpha = 3, theta = 2){
  return(rgamma(n, shape = alpha, scale = theta))
}


###### Question 2 ######

#### Vérification graphique ####
# Package de couleurs
library(viridisLite)
c_ <- plasma(5)

# On va simuler 10000 v.a. suivant f et g
n <- 10000

# Affichage de 2 graphiques
par(mfrow = c(1,2))

# Application de l'algorithme de la fonction inverse 
g <- rgen_g(n)
# histogramme de l'échantillon de v.a. suivant g
hist(g, freq = FALSE, main = "Distribution de l'échantillon", ylab = "Fréquence",xlab = "y", ylim = c(0,0.30), breaks = 40,col = "grey70", border = "grey40")
# Superposition de la courbe de la densité g
lines(x <- seq(min(g),max(g), 0.001), densite_g(x), col = c_[4], lwd = 3)
# Légende
legend("topright", "Densité g", col = c_[4], lwd = 3, box.lty = 1, inset = c(0.10,0), bg="Grey95", box.col = "Grey50")
# Application de la fonction rgamma
f <- rgen_f(n)
# Histogramme de l'échantillon de v.a. suivant f
hist(f, freq = FALSE, main = "Distribution de l'échantillon", ylab = "Fréquence",xlab = "x", ylim = c(0,0.30), breaks = 40,col = "grey70", border = "grey40")
# Superposition de la courbe de la densité f
lines(x <- seq(min(f),max(f), 0.001), densite_f(x), col = c_[3], lwd = 3)
# Légende
legend("topright", "Densité f", col = c_[3], lwd = 3, box.lty = 1, inset = c(0.10,0), bg="Grey95", box.col = "Grey50")



###### NOTA BENE ######
# Fonction qui prend une fonction h et un niveau de confiance en argument
# et renvoie un tableau contenant l'estimateur, son erreur quadratique moyenne 
# et un intervalle de confiance asymptotique bilatéral pour delta au niveau 1-alpha
# Par défaut on pose alpha = 5% pour avoir un intervalle de confiance au niveau 95%
MC <- function(h, alpha = 0.05){
  hbar <- mean(h) # estimateur 
  var <- var(h) # erreur quadratique moyenne (variance de l'estimateur) 
  IC1 <- round(mean(h)-qnorm(1-alpha/2)*sqrt(var(h)/n),8) # borne inférieure de l'IC
  IC2 <- round(mean(h)+qnorm(1-alpha/2)*sqrt(var(h)/n),8)# borne supérieure de l'IC
  # création d'un data.frame (tableau)
  d1 <- data.frame(c("", "" ,IC1)) # première colonne avec uniquement IC1
  d2 <- data.frame(c(hbar, var ,IC2)) # seconde colonne avec hbar, var et IC2
  D <- cbind(d1,d2) # fusion des colonnes
  # renommage des lignes
  rownames(D) <- c("Estimateur", "Erreur quadratique moyenne", "Intervalle de confiance")
  # renommage des colonnes
  names(D) <- c("", "")
  return(D)
}

###### Question 4 ######

#### Estimation de delta n°1 ####
p_n <- function(n, t = 0){
  X <- rgen_f(n) # On simule X suivant f
  Y <- rgen_g(n) # On simule X suivant f
  h <- 1*(X+Y>t) # On calcule h définie par l'indicatrice telle que X+Y>t
  return(h) # on renvoit la fonction h
}
#On veut une estimation pour n = 1000
n <- 1000
# On applique la fonction MC pour renvoyer le tableau regroupant l'estimateur,
# son erreur quadratique moyenne et un intervalle de confiance de delta au 
# niveau 95%
MC(p_n(n))


###### Question 5 ######

#### Vérification empirique de l'hypothèse de régime asymptotique de p_n ####

# prenons n = 1000
n <- 1000
# p_hat sera une vecteur qui sera rempli de n simulations de p_n
p_hat <- rep(0,n) # on part de p_hat un vecteur rempli de 0
for (i in 1:n){ # pour i variant de 1 à n, on modifie les n éléments de p_hat
  p_hat[i] <- mean(p_n(n)) # le i-ème élément de p_hat est modifié par p_n(n) 
}
# On applique la formule du TCL : (p_n - |E[p_1])/√\V(p_1)
TCL <- (p_hat - mean(p_hat))/sqrt(var(p_hat)) 

# Affiche d'un seul graphique 
par(mfrow = c(1,1))
# Histogramme de l'échantillon nommé TCL
hist(TCL, freq = FALSE, main = "Régime asymptotique", ylab = "Fréquence",xlab = "x", ylim = c(0,0.50), xlim = c(-3,3), breaks = 20,col = "grey70", border = "grey40")
# Superposition de la courbe de la densité d'une loi normale centrée et réduite 
lines(x<-seq(-3,3,0.01),dnorm(x),col = c_[2], lwd = 3)
# Légende
legend("topright", "Densité d'une N(0,1)", col = c_[2], lwd = 3, box.lty = 1, inset = c(0.01,0.01), bg="Grey95", box.col = "Grey50")


###### Question 7 ######
#### Estimation de delta n°2 ####

# Fonction de répartition F d'une loi gamma de paramètre de forme alpha = 3 
# et d'échelle theta = 2
fdr_F <- function(x, alpha = 3, theta = 2){
  return(pgamma(x, shape = alpha, scale = theta))
}

# Estimation de delta_n
delta_n <- function(n, t = 0){
  Y <- rgen_g(n) # On simule n v.a. suivant g
  h <- 1 - fdr_F(t - Y) # On calcule h_tilde définie par 1 - F(t - Y)
  return(h) # on renvoit la fonction h_tilde
}
#On veut une estimation pour n = 1000
n <- 1000
# On applique la fonction MC pour renvoyer le tableau regroupant l'estimateur,
# son erreur quadratique moyenne et un intervalle de confiance de delta au 
# niveau 95%
MC(delta_n(n))



###### Question 10 ######
#### Estimation de delta n°3 ####

# Estimation de delta_A_n (méthode de la variable antithétique)
delta_A_n <- function(n, t = 0){
  U <- runif(n) # on simule n v.a. uniforme 
  # on calcule h_tilde(G_fleche(U))
  h <- 1- fdr_F(t-G_fleche(U)) 
  # on calcule h_tilde(G_fleche(A(U)) avec A : x |-> 1-x
  h_A <- 1 - fdr_F(t-G_fleche(1-U)) 
  # on fait la moyenne entre h_tilde(G_fleche(U)) et h_tilde(G_fleche(A(U)))
  H <- (h + h_A)/2 
  return(H) # on renvoie (h•G_fleche(U)+h•G_fleche•A(U))/2
}
#On veut une estimation pour n = 1000
n <- 1000
# On applique la fonction MC pour renvoyer le tableau regroupant l'estimateur,
# son erreur quadratique moyenne et un intervalle de confiance de delta au 
# niveau 95%
MC(delta_A_n(n))



###### Question 11 ######

# Fonction h_tilde définie par y |--> 1-F(t-y) où F est la f.d.r. de f
h_tilde <- function(y, t = 0){
  return(1- fdr_F(t-y))
}

# Fonction h_0,1 définie par y |--> y
h_0_1 <- function(y){
  return(y)
}

# Fonction h_0,2 définie par y |--> 1_{y ≥ t - q_{epsilon}
h_0_2 <- function(y, alpha = 3, theta = 2, t = 0, epsilon = 0.6){
  q <- qgamma(epsilon, shape = alpha, scale =theta)
  return(1*(y >= t - q))
}

# Fonction permettant de calculer rho(h,h_0) (pour une v.a. Y)
rho <- function(Y,h,h_0){
  Cov_h_h_0 <- cov(h, h_0) # on calculer la covariance entre h(Y) et h_0(Y)
  var_h_0 <- var(h_0) # on calcule la variance de h_0(Y)
  var_h <- var(h) # on calcule la variance de h(Y)
  # on renvoie rho(h,h_0)
  return(Cov_h_h_0/sqrt(var_h*var_h_0))
}

# On vérifie qu'elle est le rho le plus proche de 1 entre 
# rho(h_tilde,h_0_1) et rho(h_tilde,h_0_2) 

# on va simuler n = 1000 v.a. suivant la loi de Y (densité g)
n <- 1000
# on va comparer les rho sur I = 100 itérations
I <- 100

# on va compter le nombre de fois où rho(h_tilde,h_0_1) est plus 
# grand que rho(h_tilde,h_0_2)
rho_1 <- numeric(I) # on crée un vecteur qui sera composé de I = 100 simulations de rho(h_tilde,h_0_1)
rho_2 <- numeric(I) # on crée un vecteur qui sera composé de I = 100 simulations de rho(h_tilde,h_0_2)
k <- 0 # on initialise le compteur k à 0
for (i in 1:I){ # pour i variant de 1 à I = 100,
  Y <- rgen_g(n) # on simule n v.a. suivant g
  rho_1[i] <- rho(Y,h_tilde(Y),h_0_1(Y)) # on remplit rho_1
  rho_2[i] <- rho(Y,h_tilde(Y),h_0_2(Y)) # on remplit rho_2
  # on ajoute au compteur 1 si rho(h_tilde,h_0_1) > rho(h_tilde,h_0_2) et 0 sinon
  k <- k + (rho_1[i] > rho_2[i])
}
# on renvoie le compteur k 
k
# on renvoie les valeurs moyennes de rho(h_tilde,h_0_1) et rho(h_tilde,h_0_2)
c(mean(rho_1),mean(rho_2))



# calcul du b optimal pour l'estimateur de la méthode de la variable de contrôle 
b_star <- function(Y){
  # on applique la formule donnée précédemment
  return(cov(h_tilde(Y),h_0_1(Y))/var(h_0_1(Y)))
}
Y <- rgen_g(n)
b_star(Y)


###### Question 12 ######
#### Estimation de delta n°4 ####
# Estimation de delta_cont (méthode de la variable de contrôle)
delta_cont <- function(n, mu = -13, b=2, alpha = 3, theta = 2, epsilon = 0.6, t=0 ){
  Y <- rgen_g(n) # on simule n v.a. suivant g
  b_cont <- b_star(Y) # on calcule le b optimal
  m <- mu + b/6 # on calcule m, l'espérance de Y (calculée théoriquement)
  h <- h_tilde(Y) - b_cont*(h_0_1(Y) - m) # on calcule h_cont 
  return(h)# on retourne h_cont
}

#On veut une estimation pour n = 1000
n <- 1000
# On applique la fonction MC pour renvoyer le tableau regroupant l'estimateur,
# son erreur quadratique moyenne et un intervalle de confiance de delta au 
# niveau 95%
MC(delta_cont(n))


###### Question 14 ######

#### Vérification que Sigma est inversible ####
# Valeurs des paramètres
b <- 2
mu <- -13
epsilon <- 0.6
q <- qgamma(epsilon,shape = 3, scale = 2)
t <- 0

# On calcule les éléments (variances et covariance) de Sigma à l'aide de leur formule 
# théorique qu'on a calculé
Cov_Y_Ind <- 1/2*exp(-(t-q-mu)/b)*(t-q-mu+5*b/6)
Var_Ind <- 1/2*exp(-(t-q-mu)/b)*(1-1/2*exp(-(t-q-mu)/b))
Var_Y <- 47*b*b/36

# On crée la matrice Sigma
Sigma <- matrix(c(Var_Y,Cov_Y_Ind,Cov_Y_Ind,Var_Ind),2)
Sigma

# On calcule le déterminant de Sigma
det(Sigma)

# On calcule l'inverse de Sigma
Sigma_inv <- solve(Sigma)
Sigma_inv



###### Question 17 ###### 
#### Estimation de delta n°5 ####

# Espérance de h_0_1(Y) calculée théoriquement
esp_h_0_1 <- function(mu = -13,b=2){
  return(mu + b/6)
}

# Espérance de h_0_2(Y) calculée théoriquement
esp_h_0_2 <- function(mu = -13, b=2, epsilon=0.6, t=0, alpha=3, theta=2){
  q <- qgamma(epsilon, shape = alpha, scale = theta)
  return(exp(-(t - q - mu)/b)/2)
}

# Calcul du vecteur C 
vecteur_C <- function(Y){
  Cov_1 <- cov(h_tilde(Y),h_0_1(Y)) # on calcule la covariance entre h_tilde(Y) et h_0_1(Y)  
  Cov_2 <- cov(h_tilde(Y),h_0_2(Y)) # on calcule la covariance entre h_tilde(Y) et h_0_2(Y)  
  return(c(Cov_1,Cov_2))# on renvoie le vecteur C  
}

# Estimation de delta (méthode avec beta*)
delta_beta <- function(n, Sigma){
  Y <- rgen_g(n) # on simule n v.a. suivant g
  m1 <- esp_h_0_1() # on calcule m1, l'espérance de h_0_1(Y) (calculée théoriquement)
  m2 <- esp_h_0_2() # on calcule m2, l'espérance de h_0_2(Y) (calculée théoriquement)
  # on crée le vecteur du produit scalaire avec beta*
  vect <- matrix(c(h_0_1(Y) - m1, h_0_2(Y) - m2),ncol = 2) 
  beta <- solve(Sigma)%*%vecteur_C(Y) # on calcule beta*
  h <- h_tilde(Y) - vect%*%beta # on calcule h : Y |--> 1 - F(t-Y) - < beta* , vect >
  return(h) # on renvoie h
}

#On veut une estimation pour n = 1000
n <- 1000
# On applique la fonction MC pour renvoyer le tableau regroupant l'estimateur,
# son erreur quadratique moyenne et un intervalle de confiance de delta au 
# niveau 95%
MC(delta_beta(n,Sigma))




###### Question 20 ###### 

#### Estimation de delta n°6 ####

# Estimation de delta (méthode de stratification)
Delta_strat <- function(n){
  S <- 0 # on initialise la somme à 0
  for (k in 1:n){ # Pour i variant de 1 à n
    Y <- G_fleche((k-1+runif(n))/n) # on simule n v.a. suivant G_fleche((k-1+U)/n)
    S <- S+ h_tilde(Y) # on ajoute à S h_tilde(G_fleche((k-1+U)/n))
  }
  return(S/n) # on renvoie la somme divisée par n (car Delta_n = 1/n*Somme(h_tilde(G_fleche((k-1+U)/n))))
  
}

#On veut une estimation pour n = 1000
n <- 1000
# On applique la fonction MC pour renvoyer le tableau regroupant l'estimateur,
# son erreur quadratique moyenne et un intervalle de confiance de delta au 
# niveau 95%
MC(Delta_strat(n))




###### Question 21 ###### 

#### Comparaison de la vitesse d'exécution ####
# Import de la librairie Microbenchmark
library(microbenchmark)

# On va simuler n = 100 v.a. pour chaque estimation et comparer le temps d'exécution
n <- 100
# On compare les temps de calculs
comparaison <- microbenchmark(p_n(n), delta_n(n),delta_A_n(n),delta_cont(n),delta_beta(n,Sigma),Delta_strat(n))
print(comparaison,signif=3)

#### Comparaison des variances deux à deux ####
# On va simuler n = 1000 v.a. pour chaque estimation et comparer les variances
n <- 1000
# Vecteur des variances des estimateurs
VAR <- c(var(p_n(n)), var(delta_n(n)),var(delta_A_n(n)),var(delta_cont(n)),var(delta_beta(n,Sigma)),var(Delta_strat(n)))
# on crée une matrice (ou data.frame) pour comparer les matrices deux à deux
Tb_var <- data.frame(outer(VAR,VAR,"/")) # cette matrice fait le rapport des variances deux à deux
# on renomme les lignes et les colonnes du data.frame
names(Tb_var) <- c("p_n","delta_n", "delta_A_n", "delta_cont", "delta_beta", "Delta_strat")
rownames(Tb_var) <- c("p_n","delta_n", "delta_A_n", "delta_cont", "delta_beta", "Delta_strat")

# Affichage du tableau
Tb_var


# Efficacités relatives de Delta_n par rapport aux autres estimateurs
facteurs_cout <- c(100,60,40,25,10) # facteurs de coûts par rapport aux autres estimateurs
Tb_var[6,1:5]*facteurs_cout # dernière ligne multipliée par les facteurs de coût 






####### EXERCICE 2 #######


###### Rappel des fonctions ######

# Densité f d'une loi gamma de paramètre de forme alpha = 3 et d'échelle theta = 2
densite_f <- function(x,alpha = 3,theta = 2){
  return(dgamma(x, shape = alpha, scale = theta))
}

# Densité g de paramètre mu = -13 et b = 2
densite_g <- function(y, mu = -13, b = 2){
  return((y-mu+2*b)/(4*b*b) * (y >= mu-2*b) * (y <= mu) + exp(-(y-mu)/b)/(2*b) * (y > mu))
}

# Fonction de simulation de n v.a. de densité f (via la fonction de R rgamma)
rgen_f <- function(n = 10000, alpha = 3, theta = 2){
  return(rgamma(n, shape = alpha, scale = theta))
}

# Fonction de répartition de G calculée théoriquement
fdr_G <- function(y,mu = -13, b = 2){
  return((y - (mu-2*b))^2/(8*b*b)*(mu - 2* b <= mu)*(y <= mu)  +  (1-exp(-(y-mu)/b)/2)*(y>mu))
}

# Fonction inverse généralisée de G calculée théoriquement
G_fleche <- function(u,mu = -13, b = 2){
  return((2*b*(sqrt(2*u) - 1) + mu)*(0 <= u)*(u <= 1/2)  +  (mu-b*log(2*(1-u)))*(u > 1/2)*(u < 1) )
}

# Fonction de simulation de n v.a. de densité g (via la la méthode de la fonction inverse)
rgen_g <- function(n = 10000, mu = -13, b = 2){
  # On simule n v.a. uniforme sur [0,1]
  U <- runif(n)
  # On retourne G_flèche(U) car si U ~ U([0,1]), alors G_flèche(U) et Y sont de même loi
  return(G_fleche(U))
}


###### Question 1 ######

# Densité Psi 
densite_psi <- function(x,m,k){
  return(factorial(m)/(factorial(k-1)*factorial(m-k))* 
           fdr_G(x)^(k-1)*(1-fdr_G(x))^(m-k)*densite_g(x))
}

# Probabilité d'acceptation 1/M calculée théoriquement
p_accept <- function(m, k){
  M <- factorial(m)/(factorial(k-1)*factorial(m-k)) *((k-1)/(m-1))^(k-1) *((m-k)/(m-1))^(m-k)
  return(1/M)
}

# Algorithme du rejet pour simuler n variables aléatoire suivant Psi
rejet_Psi <- function(n = 10000, m, k){
  res <- c() # vecteur résultat à remplir
  p <- p_accept(m,k) #probabilité d'acceptation calculée précédemment
  compteur <- 0 # compteur du nombre de variables aléatoires simulées
  while (length(res) < n){# tant que le vecteur résultat n'est pas de taille n,
    # on réitère le schéma suivant : 
    
    # on rajoute au compteur le nombre de v.a. simulées dans la boucle 
    compteur <- compteur+(n - length(res))%/%p +1
    
    # on simule (n - length(res)) %/% p + 1 v.a. suivant la densité g
    Y <- rgen_g((n - length(res)) %/% p +1)
    
    # on simule (n - length(res)) %/% p + 1 v.a. uniformes sur [0,1]
    U <- runif((n - length(res))%/%p +1)
    
    # On récupère uniquement les éléments du vecteur Y qui vérifie U < Psi(Y)/(M*g(Y))
    res <- c(res,Y[U < p* densite_psi(Y,m,k)/(densite_g(Y))])
  }
  # On retourne le résultat
  return(res[1:n])
  # NB : pour renvoyer la probabilité d'acceptation de l'algorithme, on return(n/compteur) 
}

###### Question 2 ######

#### Vérification graphique ####

# Package de couleurs
library(viridisLite)
c_ <- viridis(5)

# On va simuler 10000 v.a. suivant Psi
n <- 10000
# Paramètres m et k données 
m <- 15
k <- 7
# Application de l'algorithme de rejet 
x <- rejet_Psi(n,m,k)

# Histogramme de l'échantillon
hist(x, freq = FALSE,main = "Distribution de l'échantillon", ylab = "Fréquences", ylim = c(0,1), breaks = 30, col = "grey70", border = "grey40")
# Superposition de la courbe de la densité Psi
lines(y <- seq(min(x),max(x),0.001), densite_psi(y,m,k), col = c_[3], lwd = 4)
# Légende
legend("topright", "Densité Psi", col = c_[3], lwd = 4, box.lty = 1, bg = "Grey95", inset = 0.005, box.col = "grey50")






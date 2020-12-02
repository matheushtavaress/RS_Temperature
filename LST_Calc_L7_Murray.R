library(maptools)
library(raster)
library(readr)
library(rgdal)
library(rgeos)
library(stringr)

# função para calcular a proporção de vegetação no pixel, de 0 a 1 (Carlson and Ripley, 1997)
Pv = function(x) {
  Pv = ((x - NDVIs)/(NDVIv - NDVIs))^2
  return(Pv)
}

# Função para fazer a máscara de pixels
masq = function(band, pqb) {
  aux1 = values(band)
  aux2 = values(pqb)
  aux1 = aux1*aux2
  
  values(band) = aux1
  return(band)
}

# Função para fazer calcular a temperatura no trecho de rio
Tstretch = function(Tshp) {
  num = 0
  den = 0
  aux = values(Tshp)
  for (i in 1:length(aux)) {
    if (is.na(aux[i])==F) {
      num = num + aux[i]
      den = den + 1
    }
  }
  Temp = num/den
  return(Temp)
}

# Função para detectar a banda, retorna a leitura com a posição na pasta
read_band = function(list, pattern) {
  a = which(str_detect(list, pattern)==TRUE)
  band = raster(list[a])
  return(band)
}

raster_buffer = function(BQA, band, width=30) {
  prebuf = rasterToPolygons(BQA, n = 4, na.rm = TRUE, dissolve = TRUE)
  w = -width
  buf = buffer(prebuf, width = w)
  cr = crop(BQA, buf, snap = 'in') # Fazendo o crop do raster pelo shape
  fr = rasterize(buf, cr)    # Rasterizando a imagem resultante
  BQA = mask(x = cr, mask = fr)  # Incluíndo os valores no raster resultante
  BQA = extend(BQA, band)
  return(BQA)
}

mainfolder = "S:/Matheus/Archives/Doutorado/SRTR/Australia"
shape = "S:/Matheus/Archives/Doutorado/SRTR/Australia/Shapes"
landsat = "S:/Matheus/Archives/Doutorado/SRTR/Australia/L7"

# Abrindo a planilha com os dados para os algoritmos
setwd(mainfolder)
plan = read_delim("Atm_par_L8_White2.csv", col_names = T, delim = ";", escape_double = FALSE, trim_ws = TRUE)
y = plan[[1]]
m = plan[[2]]
d = plan[[3]]
w1 = plan[[4]]
w2 = plan[[8]]
tau = plan[[5]]
Lu = plan[[6]]
Ld = plan[[7]]

temps = matrix(NA, nrow = length(y), ncol = 3)
colnames(temps) = c("Tjsm", "Tjsa", "Tatc")

# Lendo o shape do trecho de rio para recorte
setwd(shape)
area = raster("Rectangle_Murray.tif")
setwd(mainfolder)
river = readOGR(dsn = "Shapes", layer = "Shape_Murray2")

# Arquivo de datas Landsat
setwd(landsat)
ListLsFiles = list.files()

for (k in seq_along(ListLsFiles)) {
  # Abrindo o diretório de cada data Landsat
  year = y[k]
  month = m[k]
  day = d[k]
  pasta = str_c(year,month,day)
  pasta = str_c(landsat,"/",pasta)
  setwd(pasta)
  
  # 1. Abrindo e tratando as imagens do Landsat para o recorte ####
  files = list.files()
  lsb2 = read_band(files, "b2.tif") # verde
  lsb6 = read_band(files, "b61.tif") # termal
  LSQ = read_band(files, "bqa.tif") # BQA
  
  # 1.1. Recorte para a área de estudo
  lsb2 = crop(lsb2, area, snap = 'in')
  lsb6 = crop(lsb6, area, snap = 'in')
  LSQ = crop(LSQ, area, snap = 'in')
  
  # 1.2. Gerando a máscara de pixels (NA = ruim, 1 = bom)
  val = values(LSQ)
  val = val==672 | val==1
  val = sapply(val, function(x) {
    if (x==T) return(1)
    else return(NA)})
  
  # 1.3. Verificando a necessidade de buffer
  if (sum(val, na.rm=T)<50) {
    next
    
  } else if (is.na(min(val))==T) {
    values(LSQ) = val
    LSQ = raster_buffer(LSQ, lsb6)
    
    # 1.4. Fazendo a multiplicação das imagens das bandas pela imagem de qualidade dos pixels
    lsb2 = masq(lsb2, LSQ)
    lsb6 = masq(lsb6, LSQ)
  }
  
  pasta = str_c(year,month,day)
    
  # 2. Carregando constantes ####
  lbdeff = 11.269 # micrômetros, valor do Lambda efetivo (comprimento de onda efetivo da banda 6 do Landsat 7)
  c1 = 1.19104*10^8 # constante de Planck
  c2 = 1.43877*10^4 # constante de Planck
  
  # Constantes termais de calibração do ETM+ (Landsat 7 User Handbook)
  K1 = 666.09 # W/m²
  K2 = 1282.71 # K
  
  # 2.1. Carregando os valores constantes da imagem do arquivo MTL
  a = which(str_detect(files, "MTL.txt")==TRUE)
  MTL = read.table(files[a], header = FALSE, quote = " ", sep = ";", dec = ".")
  a = as.character(MTL[[98,1]]) # Posição temporária do Lmax
  LMax = str_sub(a, start = 33L)  # Função para retirar apenas o valor da linha do arquivo TXT
  LMax = as.numeric(LMax) # Transformando a variável de character para numeric
  a = as.character(MTL[[99,1]]) # Posição temporária do Lmin
  LMin = str_sub(a, start = 33L)  # Função para retirar apenas o valor da linha do arquivo TXT
  LMin = as.numeric(LMin) # Transformando a variável de character para numeric
  a = as.character(MTL[[134,1]]) # Posição temporária do QcalMax
  Qmax = str_sub(a, start = 33L)  # Função para retirar apenas o valor da linha do arquivo TXT
  Qmax = as.numeric(Qmax) # Transformando a variável de character para numeric
  a = as.character(MTL[[135,1]]) # Posição temporária do QcalMin
  Qmin = str_sub(a, start = 33L)  # Função para retirar apenas o valor da linha do arquivo TXT
  Qmin = as.numeric(Qmin) # Transformando a variável de character para numeric
  
  # 2.2. Constantes para cálculo de AF1, AF2 e AF3 (Tabela de Jiménez-Muñoz et al., 2009)
  # Constantes para o cálculo apenas com w -> TIGR3 
  a = c(0.06982, -0.51041, -0.05457)
  b = c(-0.03366, -1.20026, 1.52631)
  c = c(1.04896, 0.06297, -0.32136)
  
  # 3. Transdormação de DN para temperatura de brilho ####
  # 3.1. Transformação de Digital Number (DN) em at-sensor Spectral Radiance em W/(m² sr um)
  Llbd = ((LMax - LMin)/(Qmax - Qmin))*(lsb6 - Qmin) + LMin
    
  # 3.2. Conversão para temperatura de brilho -> Chandler et al., 2009
  Tbrt = K2/log((K1/Llbd) + 1) # K
  
  # 4. Cálculos de WST ####
  # 4.1. SC com w do MOD07
  if (is.na(w1[k])) {
    temps[1,k] = NA
  } else {
    AF = vector(mode = 'numeric', length = 3L)
    for (i in 1:3) {
      AF[i] = a[i]*w1[k]^2 + b[i]*w1[k] + c[i]
    }
    
    # Cálculo de gama
    gama = (Tbrt^2)/(c2*Llbd*((lbdeff^4)*Llbd/c1 + (1/lbdeff)))
    
    # Cálculo de delta
    delta = Tbrt - gama*Llbd
    
    # Cálculo da temperatura de superfície para toda área
    LST = gama*(((AF[1]*Llbd + AF[2])/Emiss) + AF[3]) + delta - 273.15
    
    # Gravando a banda de temperatura
    name = str_c("LST_SC_MOD07_", pasta, ".tif")
    writeRaster(LST, filename = name, format = "GTiff") # Salvando o raster
    
    # Recortando para o trecho de rio e escrevendo
    cr = crop(LST, river, snap = 'out') # Fazendo o crop do raster pelo shape
    fr = rasterize(river, cr)    # Rasterizando a imagem resultante
    LST = mask(x = cr, mask = fr)  # Incluíndo os valores no raster resultante
    
    temps[1,k] = Tstretch(LST)
  }
  
  # 4.2. SC com w do AtcCorr
  if (is.na(w2[k])) {
    temps[2,k] = NA
  } else {
    AF = vector(mode = 'numeric', length = 3L)
    for (i in 1:3) {
      AF[i] = a[i]*w2[k]^2 + b[i]*w2[k] + c[i]
    }
    
    # Cálculo de gama
    gama = (Tbrt^2)/(c2*Llbd*((lbdeff^4)*Llbd/c1 + (1/lbdeff)))
    
    # Cálculo de delta
    delta = Tbrt - gama*Llbd
    
    # Cálculo da temperatura de superfície para toda área
    LST = gama*(((AF[1]*Llbd + AF[2])/Emiss) + AF[3]) + delta - 273.15
    
    # Gravando a banda de temperatura
    name = str_c("LST_SC_ACPC_", pasta, ".tif")
    writeRaster(LST, filename = name, format = "GTiff") # Salvando o raster
    
    # Recortando para o trecho de rio e escrevendo
    cr = crop(LST, river, snap = 'out') # Fazendo o crop do raster pelo shape
    fr = rasterize(river, cr)    # Rasterizando a imagem resultante
    LST = mask(x = cr, mask = fr)  # Incluíndo os valores no raster resultante
    
    temps[2,k] = Tstretch(LST)
  }
    
  # 4.3. cálculo: variáveis atmosféricas do AtcCorr
  # Correção de Lsensor
  Lcorr = (Llbd - Lu[k] - tau[k]*Ld[k]*(1 - Emiss))/(Emiss*tau[k])
  
  # Cálculo da temperatura de brilho corrigida
  LST = K2/log((K1/Lcorr) + 1) - 273.15
  
  # Gravando banda de temperatura
  name = str_c("LST_ACPC_", pasta, ".tif")
  writeRaster(LST, filename = name, format = "GTiff") # Salvando o raster
  
  # 7. Recortando para o trecho de rio e escrevendo
  cr = crop(LST, river, snap = 'out') # Fazendo o crop do raster pelo shape
  fr = rasterize(river, cr)    # Rasterizando a imagem resultante
  LST = mask(x = cr, mask = fr)  # Incluíndo os valores no raster resultante
  
  temps[3,k] = Tstretch(LST)
}

# Escrevendo os resultados
setwd(mainfolder)
delta = data.frame(year, month, day, temps)
write_delim(delta, "Results_L7_Murray.csv", delim = ";", col_names = T)

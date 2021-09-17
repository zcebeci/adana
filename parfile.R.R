#Genel parametreler
#maxiter=100				#[1,Inf]			Muhtelif fonksiyonlar : Maksimum iterasyon sayısı
#g=1					#[1,maxiter], dinamik		Muhtelif fonksiyonlar : İterasyon/Kuşak no	
#gmax					# maxiter: dinamik		Muhtelif fonksiyonlar : Maksimum kuşak sayısı	
#monitorfunc				#{monprogress, montestfuncton}	adana : İzleme fonksiyonu adı			
#verbose=FALSE				#{TRUE, FALSE}			adana : Çıktı görüntüleme tercihi			
#parfile=NULL				#"gaparameters.R"		adana : Parametre dosyası (bu dosya)

# Başlatma parametreleri
#n=20					#[2,Inf]			initbin, initval, initperm : Popülasyon büyüklüğü	
#m=8					#[2,Inf]			initbin, initval, initperm : Kromozom uzunluğu	
#prevpop=NULL 				#n*m atris veya data frame	initbin, initval, initperm : Ön popülasyon
#type=1					#{1,2}				initbin, initval, initperm : Pop. matrisi türü	
#nmode="real"				#{"real", "integer"}		initval	: Değer türü			
#lb=c()					#vektör				initval	: Alt sınır değerleri vektörü			
#ub=c()					#vektör				initval	: Üst sınır değerleri vektörü			
#permset=0:9				#vektör				initperm : Permütasyon değerieri seti
	
#Seleksiyon parametreleri
selb=0.5				#[0,Inf]			selescale : Beta parametresi
selbc=0.5				#[0, 1]				selers : Taban parametresi	
selc=2					#[1, 2]				selsscale, selsscale2 :	Ölçeklendirme parametresi	
selk=1.005				#[0.95, 1.1]			selpscale : Güç parametresi			
sells=1.5				#[1.2, 2]			sellscale : ölçeklendirme parametresi			
selns=0.5				#[0, Inf]			selnlrs : Polinomiyal katsayı			
selps=0.5				#[0, 1]				seltrunc : Kesme eşiği %			
sels=1.5				#[1, 2]				sellrs,sellrs3, selrscale2 : Seçim basıncı	
selt=2					#[2, n]				seltour, seltour2: Turnuva büyüklüğü		
selt0=50				#[5, 100]			selboltour: Başlangıç sıcaklığı			
selw=2					#[2, 10]			selwscale : Pencere büyüklüğü			
reptype=FALSE				#{FALSE, TRUE}			seltour, seltour2: İadeli-İadesiz örnekleme		
#fmin=0					#[0,Inf]			selwscale: Minimum uyum, dinamik			
#selg					#[1, Inf]			selboltour: Kuşak no, dinamik	
#selgmax				#[2, Inf]			selboltour: Kuşak no, dinamik	

#Çaprazlama parametreleri
cxpc=0.9				#[0,1]				Tüm çaprazlama fonks.: Çaprazlama oranı
cxon=2	 				#[1,Inf]			Tüm çaprazlama fonks. : Çiftleşme başına yavru sayısı		
cxk=2					#[1,Inf]			kpx : Kesme noktası sayısı
cxps=0.5				#[0, 1]				hux, ux, ux2, dc: Olasılıksal eşik
cxa = 0					#[0, 1]	 			lap: Alt sınır değerleri
cxb = 0.15				#>0, 0.15, 0.35	 		lap: Üst sınır değerleri
cxealfa = 1				#[1, Inf]			elx: Genişletme oranı
#cxalfa 	    			#[0, 1]				sax, wax, ebx: Rastlantısal alfa değeri,dinamik
#cxoxk					#[1,m/2]			ox2: Rastlantısal k değeri,dinamik

# Mutasyon parametreleri
mutpm = 0.05				#[0,1]				Tüm mutasyon fonk.: Mutasyon oranı
mutb = 0.5				#(0,Inf)			nunimut, nunimut2 # Tek düzeliği önleme için üs 
mutpow = 2				#(0,Inf)			powmut, powmut2: Üs değeri
#mutmy 					#[-Inf, Inf] 			gaussmut3 : Değişken ortalamaları vektörü
#mutsdy 				#[-Inf, Inf]			gaussmut3 : Değişken std. sapmaları vektörü
#mutq 					#[0,1]				bsearchmut1 		
#lb					#[-Inf, Inf]
#ub					#[-Inf, Inf]

# Popülasyon yenileme parametreleri
reps=2					#[1,n-1]			elitism: Seçkinci yenileme birey sayısı
repk = 10				#[2,n-1]			grrobin: Karşılaştırılacak birey sayısı
reppars=c(1:2)				#c(1:n-1)			ssrfamtour, ssrx: Karşılaştırmaya alınacak birey sayısı
#lambda 				#[1,n-1]			grmuplambda2, grmuplambda3, grmuplambda4: Yavru sayısı	

#Sonlandırma parametreleri
abdif					#[-Inf, Inf]			terminate : Ortalama ve en iyi uyum arası minimum fark	
bestdif					#[0, Inf]			terminate : Son iki en iyi uyum arasındaki fark yüzdesi	
optdif=1e-06				#[0,Inf]			terminate : Global optimum değerine yaklaşma değeri
rmcnt=10				#[1,maxiter]			terminate : Son uyum ortalamaları sayısı 		
rmdif=1e-06				#[0,Inf]			terminate : on k en iyi uyum ortalamaları arası minimum fark	
phidif=1e-03				#[0,Inf]			terminate : Phi yakınsaması		
meandif=1e-06				#[0,Inf]			terminate : Son iki ort. arası minimum fark		
sddif=1e-03				#[0,Inf]			terminate : Son iki standart sapma arasındaki minimum fark		
rangedif=1e-03				#[0,Inf]			terminate : Minimum ve maksimum farkı (Değişim aralığı)		
mincv=0.05				#[0,1]				terminate : Minimum varyasyon katsayısı		
simlev=0.95				#[0,1]				terminate : Uyum değerleri benzerlik yüzdesi		
maxtime=60				#[0,Inf]			terminate : Çalışma süresi, max(60, maxiter/0.0003)
#maxiter				#[1,Inf]			terminate : Maksimum iterasyon, dinamik
#objval					#[-Inf, Inf]			terminate : Global optimum değeri

#Melez GA parametreleri
hgastep=10				#[1,maxiter]			hgafunc: Melezleme adım sayısı 		
hgans=2					#[1,n]				hgafunc: Melezlenecek birey sayısı	
hgaftype="r"				#{"r", "w", "b"}		hgafunc: Uyum değeri türü

#Uyarlamalı GA parametreleri
cxpc2					#[0,1], dinamik			leitingzhi: Uyarlanmış çaprazlama oranı
mutpm2					#[0,1], dinamik			leitingzhi: Uyarlanmış mutasyon oranı	
adapa=0.7				#[0,1]				leitingzhi: Uyarlama eşiği 1
adapb=0.5				#[0,1]				leitingzhi: Uyarlama eşiği 2
adapc=0.2				#[0,1]				adana3: Uyarlama eşiği 1
adapd=0.05				#[0,1]				adana3: Uyarlama eşiği 2


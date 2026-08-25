program edson
    implicit none

integer, parameter :: L = 512
integer, parameter :: N = L*L
integer, parameter :: t = 50000000
integer, parameter :: NRK = 25
integer, parameter :: NE =100*L**2
double precision, parameter :: EXPMUT = 4.0
double precision :: H = 0.01

integer :: STA, VZS, ESP, GENE
double precision :: PMF=0.01
double precision :: PMV, PRM, PRA, NORM, PROBMUT_A
double precision :: NRKK = NRK
double precision :: MUT = 10**(-EXPMUT)

integer :: rede(L,L), especie(N), unic_ESP, tamanho_especie(NE)
double precision :: recurso(N,NRK), KMS(NE,NRK), PR(NRK)




! Inter�veis 
integer :: i, j, k

! Contadores
integer :: cont_t, cont_v ! contador de tempo e de visitas
double precision :: cont_1, cont_2


! ARQUIVOS 
integer :: v1, hv1
character*200 entrada, Saida




!------------------------------ Inicializa��o dos par�metros  --------------------------------
!Abaixo � inicializada a matriz 'KMS' com valores aleat�rios para a esp�cie inicial (ESP).
ESP=1
do i = 1, NRK
    KMS(ESP,i)= Box_Muller()
!     write(*,*) KMS(ESP,i)
end do

!Aqui, � calculada a soma dos valores na primeira linha da matriz 'KMS' e armazenada em 'cont_1'.
cont_1 = 0
do i = 1, NRK
    cont_1 = cont_1 + KMS(1, i)
!    write(*,*) cont_1 
end do

!Aqui � verificado se a soma dos valores da primeira linha de 'KMS' � diferente do n�mero de recursos 'NRKK'. Se for diferente, calcula o fator de normaliza��o 'NORM'.
if (cont_1 /= NRKK) then
    NORM = NRKK / cont_1
!    write(*,*) NORM 
end if

!Essas linhas normalizam os valores na primeira linha da matriz 'KMS' multiplicando pelo fator 'NORM' e recalculando 'cont_2' como a soma desses valores normalizados.
cont_2 = 0
do i = 1, NRK
    KMS(1, i) = KMS(1, i) * NORM
    cont_2 = cont_2 + KMS(1, i)
!   write(*,*) cont_2, KMS(1, i)
end do


!------------------------------------------ ARQUIVOS ----------------------------------------!

!Parte do programa que grava os relevos em um arquivo e salva posteriormente, abrindo e fechando os arquivos.
do i = 1, NRK
    100 format (a, i6, a, f8.2, a, i6, a)
    write (entrada, 100) "Relevo-L=", L, "-H=", H, "-", i, ".dat"
    
    Saida = " "
    hv1 = 0
    do v1 = 1, len(entrada)
        if (entrada(v1:v1) /= " ") then
            hv1 = hv1 + 1
            Saida(hv1:hv1) = entrada(v1:v1)
        end if
    end do
    
    open (400, file = Saida, status='old')
    do j = 1, N
        read(400,*) recurso(j, i)
!        write(*,*)i,j,recurso(j,i)
    end do
    close(400)
end do
!___________________________________Abrir arquivos para escrita___________________________________________________________
open(unit= 1, file = "N_ESP.dat", status='unknown')
open(unit= 2, file = "KMS_ESP_SOBRE.dat", status='unknown')
open(unit= 3, file = "ESP_SIT.dat", status='unknown')


!______________________________________________________________________________________________


do i = 1, N
    especie(i) = ESP
!    write(*,*) especie(i), i
end do

cont_t = 0
cont_v = 0

! -------------------------------------------------------------------------------------------------------!
!                                         Aplica��o da Din�mica no Modelo                                !
! -------------------------------------------------------------------------------------------------------!


1 continue  ! Ponto de in�cio do ciclo principal XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!Sele��o aleat�ria de um s�tio na rede.
STA = AleaN(N)
!write(*,*) STA

!Este trecho verifica se a esp�cie no s�tio 'STA' � igual a zero. Se for, incrementa 'cont_v' em 1. 
if (especie(STA) == 0) then
    cont_v = cont_v + 1
!    write(*,*) cont_v, cont_t
    goto 2  ! Se a esp�cie no local STA � zero, vai para o ponto 2 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
end if

!Aqui, 'cont_v' est� sendo incrementado em 1,  para contar o n�mero de vezes que este bloco de c�digo foi executado.
cont_v = cont_v + 1

!A vari�vel 'PMV' est� sendo atribu�da a um valor aleat�rio usando a fun��o 'Alea()'. � usado para simular uma probabilidade de morte vari�vel.
PMV = Alea()

!write(*,*) PMV


!Este trecho verifica se 'PMV' � menor ou igual a 'PMF' (probabilidade de morte fixa). Se isso for verdadeiro, a esp�cie no s�tio 'STA' � definida como zero, simbolizando a morte dessa esp�cie no s�tio. 
if (PMV <= PMF) then
    especie(STA) = 0
!    write(*,*) PMV, PMF, especie(STA)
    goto 2  ! Se a esp�cie morre, volta ao CONTADOR DE VISITAS   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
end if


!Este trecho faz parte de um loop que verifica se h� alguma esp�cie vizinha ao s�tio 'STA' que seja igual a zero. Se encontrar um vizinho com esp�cie zero, o c�digo vai para o ponto "3" no c�digo. Caso contr�rio, ele continua no loop.
do i = 1, 4
    !write(*,*) especie(viz(L, STA, i)), i
    if (especie(viz(L, STA, i)) == 0) then
        goto 3  ! Se n�o tem esp�cie, vai para o ponto 3 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    end if
end do

goto 2  ! Se todos os vizinhos est�o ocupados, volta ao CONTADOR DE VISITAS XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

3 continue  ! Ponto de in�cio da reprodu��o XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!write(*,*) "PASSOU"
!STOP

!Aqui � calculada a probabilidade de reprodu��o para cada recurso do 'STA'.
do i = 1, NRK
    PR(i) = recurso(STA, i) / (KMS(especie(STA), i) + recurso(STA, i))
!    write(*,*) PR(i), KMS(especie(STA), i), recurso(STA, i)
end do
!stop

!Aqui, 'PRM' recebe o valor m�nimo das probabilidades calculadas anteriormente no vetor 'PR'. Ou seja, 'PRM' representa a menor probabilidade entre os recursos dispon�veis para a esp�cie no s�tio 'STA'.
PRM = minval(PR)
!write(*,*) PRM
!PAUSE
!stop

!'PRA' � atribu�do a um valor aleat�rio usando a fun��o 'Alea()' para simular uma probabilidade de reprodu��o aleat�ria.
PRA = Alea()
!write(*,*) PRA
!PAUSE


!Verificamos se a probabilidade aleat�ria 'PRA' � maior que a menor probabilidade 'PRM' calculada anteriormente. Se essa condi��o for verdadeira, isso significa que a esp�cie no s�tio 'STA' n�o ir� reproduzir.
if (PRA > PRM) then
!write(*,*) PRA, PRM
!PAUSE
    goto 2  ! Se n�o reproduziu, vai para o ponto 2 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
end if



4 continue  ! Ponto de in�cio da busca por vizinho vazio XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! 'VZS' � a escolha de um vizinho aleat�rio entre 1 e 4.
VZS = AleaN(4)

!write(*,*) VZS
!PAUSE


!Este trecho verifica se a esp�cie no vizinho selecionado ('viz(L,STA,VZS)') � diferente de zero. Se a condi��o for verdadeira, ou seja, se houver uma esp�cie no vizinho, o programa imprime "TEM" (comentado) e retorna ao ponto "5" no c�digo. Isso significa que o programa continua a selecionar aleatoriamente um vizinho at� encontrar um vizinho sem esp�cie.
if (especie(viz(L, STA, VZS)) /= 0) then
!write(*,*) especie(viz(L, STA, VZS))
!PAUSE
    goto 4  ! Se o vizinho n�o est� vazio, continua buscando XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
end if


!A vari�vel 'PROBMUT_A' � atribu�da a um valor aleat�rio usando a fun��o 'Alea()'. Isso � usado para simular uma probabilidade aleat�ria.
PROBMUT_A = Alea()

!write(*,*) PROBMUT_A
!PAUSE


!Esta estrutura verifica se o valor de muta��o 'MUT' � menor que a probabilidade aleat�ria 'PROBMUT_A'. Se isso for verdadeiro, significa que ocorrer� uma muta��o na esp�cie. O c�digo ent�o atribui a esp�cie no vizinho 'viz(L,STA,VZS)' � esp�cie no s�tio 'STA', simulando assim uma muta��o.
if (MUT < PROBMUT_A) then
    especie(viz(L, STA, VZS)) = especie(STA)
    goto 2  ! Se N�O houve muta��o, volta ao CONTADOR DE VISITAS  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
end if

ESP = ESP + 1

!Nesta linha, a esp�cie no vizinho 'viz(L,STA,VZS)' � alterada para a esp�cie atual 'ESP' mais um. Isso sugere um processo de evolu��o ou mudan�a de esp�cies em um determinado s�tio.
especie(viz(L, STA, VZS)) = ESP 
!write(*,*) especie(viz(L, STA, VZS))
!PAUSE



!Neste loop, as constantes 'KMS' para a nova esp�cie no vizinho 'viz(L,STA,VZS)' s�o definidas iguais �s constantes da esp�cie atual 'STA'. Isso indica que as constantes de muta��o ou adapta��o da nova esp�cie s�o as mesmas que as da esp�cie atual.
do i = 1, NRK
    KMS(especie(viz(L, STA, VZS)), i) = KMS(especie(STA), i)
!    write(*,*) KMS(especie(viz(L, STA, VZS)), i)
end do
!PAUSE

!Aqui, 'GENE' � atribu�do a um valor aleat�rio entre 1 e 'NRK' usando a fun��o 'AleaN()'. Significa uma sele��o aleat�ria de um gene dentro das caracter�sticas da esp�cie no vizinho 'viz(L,STA,VZS)'.
GENE = AleaN(NRK)

!write(*,*) GENE
!PAUSE

!write(*,*) KMS(especie(viz(L, STA, VZS)), GENE)

!write(*,*) "--------"

!Nesta linha, a fun��o 'Box_Muller2()' � usada para aplicar a transforma��o ao valor da constante 'KMS' do gene selecionado ('GENE') na esp�cie do vizinho 'viz(L,STA,VZS)'.
KMS(especie(viz(L, STA, VZS)), GENE) = Box_Muller2(KMS(especie(viz(L, STA, VZS)), GENE))

!write(*,*) KMS(especie(viz(L, STA, VZS)), GENE)
!PAUSE

!Este loop imprime as constantes 'KMS' da esp�cie no vizinho 'viz(L,STA,VZS)'. Pode ser usado para verificar as mudan�as nas constantes ap�s a aplica��o da transforma��o.
do i = 1, NRK
!     print*, KMS(especie(viz(L, STA, VZS)), i)
end do
!stop
!PAUSE

cont_1 = 0

!Este loop soma os valores das Constantes de Meia e Satura��o ('KMS') de uma esp�cie no vizinho 'viz(L,STA,VZS)'. 'especie(viz(L,STA,VZS))' representa a esp�cie nesse vizinho e 'i' itera sobre os recursos ('NRK'). Portanto, 'cont_1' � a soma desses valores para todos os recursos.
do i = 1, NRK
    cont_1 = cont_1 + KMS(especie(viz(L, STA, VZS)), i)
!    write(*,*)cont_1
end do
!stop
!Aqui, 'NORM' � calculado como o n�mero total de recursos 'NRKK' dividido pela soma dos valores de 'KMS' da esp�cie no vizinho. Trata-se de uma normaliza��o onde o n�mero de recursos dispon�veis � ajustado com base na soma das 'KMS' dessa esp�cie no vizinho.
NORM = NRKK / cont_1

!write(*,*) NORM
!PAUSE

!write(*,*) "------"

!Este loop calcula as novas ('KMS') para a esp�cie no vizinho 'viz(L,STA,VZS)' multiplicando os valores antigos de 'KMS' pelo fator de normaliza��o 'NORM'. 
do i = 1, NRK
    KMS(especie(viz(L, STA, VZS)), i) = KMS(especie(viz(L, STA, VZS)), i) * NORM
!write(*,*) KMS(especie(viz(L, STA, VZS)), i) 
end do


cont_2 = 0
!Aqui, 'cont_2' � inicializado como zero e, em seguida, � somado com os novos valores de 'KMS' da esp�cie no vizinho 'viz(L,STA,VZS)'. 
do i = 1, NRK
    cont_2 = cont_2 + KMS(especie(viz(L, STA, VZS)), i)
!    write(*,*)cont_2
end do
!stop


2 continue  ! Ponto de verifica��o da contagem de visitas XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX




if (cont_v == N) then
    
    
cont_t = cont_t + 1
cont_v = 0

if ( mod(cont_t,250000)==0 ) then
    unic_ESP = unic(especie)
write(1,*) cont_t, unic_ESP
!write(*,*)cont_t,unic_ESP
    
end if
! Contabilizar as esp�cies por tempo------------------------------------------------------------



!--------------------------------------------------------------------------------------------- 
if (cont_t == t) then
   
! Imprimir as esp�cies, relevo e Constante de mei satura��o de cada relevo  
do i = 1, ESP
    if (any(especie == i)) then
        do j = 1, NRK
            write(2,*) i, j, KMS(i,j)
        end do
    end if
end do


! Imprimir as esp�cies em todos os s�tios
!write(3,*) "Esp�cies em cada s�tio ap�s a execu��o:"

do i = 1, N
    write(3,*) i, especie(i)
end do
    
    
    
    
    
    

!---------------------------------------------------------------------------------------------   



    
    
! ---------------------------------Fechamento dos arquivos----------------------------------------
    
close(1)  
close(2)
close(3)

stop
end if
end if



goto 1  ! Volta ao in�cio do ciclo principal XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX









! -------------------------------------------------------------------------------------------------------!
!                                         FUN��ES UTILIZADAS                                             !
! -------------------------------------------------------------------------------------------------------!



contains


!----------------Fun��o de n�mero aleat�rio de tamanho N----------------------!
function AleaN(N)
    double precision x1
    integer AleaN, N
    call random_number(x1)
    AleaN = int(x1 * N) + 1
end function AleaN


!----------------Fun��o de n�mero aleat�rio----------------------!

function Alea()
    double precision Alea, x1
    call random_number(x1)
    Alea = x1
end function Alea

 !----------------Fun��o Gaussiana---------------------!
function Box_Muller()
    double precision Box_Muller
    double precision :: x1, x2, y1, y2
    double precision :: var, pi, medfi
    var = 0.001
    pi = 3.1415
    medfi = 1
    call random_number(x1)
    call random_number(x2)
    y1 = sqrt(-2 * log(x1) * var) * cos(2 * pi * x2) + medfi
    y2 = sqrt(-2 * log(x1) * var) * sin(2 * pi * x2) + medfi
    Box_Muller = y2
end function Box_Muller

function Box_Muller2(kms2)
    double precision Box_Muller2, kms2
    double precision :: x1, x2, y1, y2
    double precision :: var, pi, medfi
    var = 0.001
    pi = 3.1415
    medfi = kms2
    call random_number(x1)
    call random_number(x2)
    y1 = sqrt(-2 * log(x1) * var) * cos(2 * pi * x2) + medfi
    y2 = sqrt(-2 * log(x1) * var) * sin(2 * pi * x2) + medfi
    if (y1 <= 0) then
        Box_Muller2 = y2
    else
        Box_Muller2 = y1
    end if
end function Box_Muller2

!----------------Fun��o para verificar as esp�cies---------------------!
function unic(vetor) result(contagem)
    integer, intent(in) :: vetor(:)
    integer :: contagem
    integer :: i, j
    logical :: valor_duplicado
    contagem = 1
    do i = 2, size(vetor)
        valor_duplicado = .false.
        do j = 1, i-1
            if (vetor(i) == vetor(j)) then
                valor_duplicado = .true.
                exit
            end if
        end do
        if (.not. valor_duplicado) then
            contagem = contagem + 1
        end if
    end do
end function unic


!----------------Fun��o VIZINHAN�A---------------------!

function viz(L, sitio, direcao) result(num_vizinho)
    integer, intent(in) :: L, sitio, direcao
    integer :: num_vizinho
    integer :: i, j
    i = (sitio - 1) / L + 1
    j = modulo(sitio - 1, L) + 1
    select case (direcao)
        case (1)
            if (i == 1) then
                num_vizinho = (L - 1) * L + j
            else
                num_vizinho = sitio - L
            end if
        case (2)
            if (j == L) then
                num_vizinho = sitio - (L - 1)
            else
                num_vizinho = sitio + 1
            end if
        case (3)
            if (i == L) then
                num_vizinho = j
            else
                num_vizinho = sitio + L
            end if
        case (4)
            if (j == 1) then
                num_vizinho = L * i
            else
                num_vizinho = sitio - 1
            end if
    end select
end function viz




end program edson

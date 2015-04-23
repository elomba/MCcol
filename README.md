

<center><h1> MCcol Project</h></center>

<p><b> Team: MC-cool </b></p>
<p> Members: </p>
<ul>
    <li>    Noé G. Almarza (noe@iqfr.csic.es)</li>
    <li>    Enrique Lomba (enrique.lomba@csic.es)</li>
    <li>    Eva G. Noya (eva.noya@iqfr.csic.es)</li>
    <li>	  Felipe J. Blas (felipe@uhu.es)</li>
</ul>
        
This repository has been created to control file versions of
 gpMC (<font color="red">a general purpose Monte Carlo</font> program).

Requisites
==========

- ifort version 13.0 or higher
- ecliple Luna

Install
=======

<p> To Download the application use </p>

<font color="blue">git clone https://github.com/elomba/MCcol.git</font> or import the URI from Eclipse

Buenas prácticas
================

<ul>
<li>Usar los estándares de F90, sobre todo para magnitudes vectoriales                                                                                                                                                           
 y módulos (para compartir datos y funciones entre distintas partes del
 programa)

<li> <font color="red"> Muy importante :</font> Al introducir cambios enmarcarlos entre comentarios comenzando por !!!InicialesVersion, eg.

<font color="blue">

!!EL0.1

    do i =1, n
    
    ....
    
    enddo
    
!!EL0.1
</font>

<li>Subir las versiones definitivas (push) a github con versiones consecutivas. Cambios menores se indican con letras (e.g. 0.0.3a). Se considera que la 1.0 será la primera versión distribuible 

<li> Las definiciones generales de variables compartidas entre los
distintos elementos (ficheros) del programa van en
Definitions.f90. Conviene modificar lo menos posible este fichero, y
cuando se haga marcar bien con comentarios la modificación y el autor

<li> Al importar un modulo a una subrutina, utilizar la directiva "only"
para importar sólo las variables necesarias y evitar potenciales
conflictos: e.g. use rundata, only : drmax 

<li> Al incorporar nuevas características al programa hay que modificar 
   necesariamente <font color="blue">Main.f90 e Init.f90</font> (este último contiene las
   subrutinas de inicialización y lectura de datos). Mi sugerencia es
   que la mayor parte de los cambios que se hagan vayan en ficheros
   separados (incluido por ejemplo nuevos procedimientos para
   generar/leer la configuración inicial), añadiendo los
   correspondientes objetos (f<font color="blue">ile.f90 --> file.o</font> en la entrada
   <font color="blue">OBJ_CORE</font> del <font color="blue">Makefile</font>. 

<li> Al definir la precisión de las variables, utilizad siempre <font color="blue">(wp)</font> del
modulo set_precision, eso nos permitirá modificar la precisión de
forma general si se quiere obtener una versión en precisión
sencilla. En principio para GPUs habrá que utilizar precisión mixta,
con los acumuladores en doble precisión.

<li> Cuando se crean subrutinas y funciones que tengan argumentos, hay
que incluir una declaración en el modulo <font color="blue">Interfaces</font> en
 <font color="blue">Definitions.f90</font>. En el subprograma en el que se haga uso de las
subrutinas/funiones hay que incluir
<font color="blue">

use interfaces, only : subroutine_name</font>

de esa forma el compilador chequea errores en el número y tipo de
argumentos usados.

</ul>

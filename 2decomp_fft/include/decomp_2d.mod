	  TV  »   k820309              11.0        ÁÚíY                                                                                                           
       decomp_2d.f90 DECOMP_2D       &       DECOMP_2D_INIT DECOMP_2D_FINALIZE DECOMP_INFO_INIT DECOMP_INFO_FINALIZE PARTITION DECOMP_2D_ABORT GET_DECOMP_INFO MYTYPE REAL_TYPE COMPLEX_TYPE MYTYPE_BYTES NX_GLOBAL NY_GLOBAL NZ_GLOBAL COMMLOC NRANK NPROC DECOMP_2D_COMM_CART_X DECOMP_2D_COMM_CART_Y DECOMP_2D_COMM_CART_Z DECOMP_INFO XSTART XEND XSIZE YSTART YEND YSIZE ZSTART ZEND ZSIZE gen@TRANSPOSE_X_TO_Y gen@TRANSPOSE_Y_TO_Z gen@TRANSPOSE_Z_TO_Y gen@TRANSPOSE_Y_TO_X gen@UPDATE_HALO gen@ALLOC_X gen@ALLOC_Y gen@ALLOC_Z                                                    
                                                             u #TRANSPOSE_X_TO_Y_REAL    #TRANSPOSE_X_TO_Y_COMPLEX 	   #         @     @X                                             #TRANSPOSE_X_TO_Y_REAL%SIZE    #TRANSPOSE_X_TO_Y_REAL%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                                 SIZE               @                                 PRESENT           
 @@                                                
 4             &                   &                   &                                                     D@@                                                
 5              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO    #         @     @X                            	                  #TRANSPOSE_X_TO_Y_COMPLEX%SIZE 
   #TRANSPOSE_X_TO_Y_COMPLEX%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                            
     SIZE               @                                 PRESENT           
 @@                                                 6             &                   &                   &                                                     D@@                                                 7              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO                                                          u #TRANSPOSE_Y_TO_Z_REAL    #TRANSPOSE_Y_TO_Z_COMPLEX    #         @     @X                                             #TRANSPOSE_Y_TO_Z_REAL%SIZE    #TRANSPOSE_Y_TO_Z_REAL%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                                 SIZE               @                                 PRESENT           
 @@                                                
 D             &                   &                   &                                                     D@@                                                
 E              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO    #         @     @X                                              #TRANSPOSE_Y_TO_Z_COMPLEX%SIZE    #TRANSPOSE_Y_TO_Z_COMPLEX%PRESENT    #SRC    #DST    #OPT_DECOMP                  @                                 SIZE               @                                 PRESENT           
 @@                                                 F             &                   &                   &                                                     D@@                                                 G              &                   &                   &                                                     
 @@                                   è             #DECOMP_INFO                                                          u #TRANSPOSE_Z_TO_Y_REAL    #TRANSPOSE_Z_TO_Y_COMPLEX !   #         @     @X                                             #TRANSPOSE_Z_TO_Y_REAL%SIZE    #TRANSPOSE_Z_TO_Y_REAL%PRESENT    #SRC    #DST    #OPT_DECOMP                   @                                 SIZE               @                                 PRESENT           
@@@                                                
 T             &                   &                   &                                                     D@@                                                
 U              &                   &                   &                                                     
 @@                                    è             #DECOMP_INFO    #         @     @X                            !                  #TRANSPOSE_Z_TO_Y_COMPLEX%SIZE "   #TRANSPOSE_Z_TO_Y_COMPLEX%PRESENT #   #SRC $   #DST %   #OPT_DECOMP &                 @                            "     SIZE               @                            #     PRESENT           
@@@                             $                    V             &                   &                   &                                                     D@@                             %                    W              &                   &                   &                                                     
 @@                              &     è             #DECOMP_INFO                                                          u #TRANSPOSE_Y_TO_X_REAL '   #TRANSPOSE_Y_TO_X_COMPLEX -   #         @     @X                           '                  #TRANSPOSE_Y_TO_X_REAL%SIZE (   #TRANSPOSE_Y_TO_X_REAL%PRESENT )   #SRC *   #DST +   #OPT_DECOMP ,                 @                            (     SIZE               @                            )     PRESENT           
 @@                             *                   
 d             &                   &                   &                                                     D@@                             +                   
 e              &                   &                   &                                                     
 @@                              ,     è             #DECOMP_INFO    #         @     @X                            -                  #TRANSPOSE_Y_TO_X_COMPLEX%SIZE .   #TRANSPOSE_Y_TO_X_COMPLEX%PRESENT /   #SRC 0   #DST 1   #OPT_DECOMP 2                 @                            .     SIZE               @                            /     PRESENT           
 @@                             0                    f             &                   &                   &                                                     D@@                             1                    g              &                   &                   &                                                     
 @@                              2     è             #DECOMP_INFO                                                           u #UPDATE_HALO_REAL 3   #UPDATE_HALO_COMPLEX ;   #         @     @X                            3                  #UPDATE_HALO_REAL%SIZE 4   #UPDATE_HALO_REAL%PRESENT 5   #IN 6   #OUT 7   #LEVEL 8   #OPT_DECOMP 9   #OPT_GLOBAL :                 @                            4     SIZE               @                            5     PRESENT           
 @@                             6                   
 z             &                   &                   &                                                   D @@                             7                   
 {              &                   &                   &                                                     
   @                              8                      @@                              9     è              #DECOMP_INFO               @@                              :            #         @     @X                            ;                  #UPDATE_HALO_COMPLEX%SIZE <   #UPDATE_HALO_COMPLEX%PRESENT =   #IN >   #OUT ?   #LEVEL @   #OPT_DECOMP A   #OPT_GLOBAL B                 @                            <     SIZE               @                            =     PRESENT           
 @@                             >                    ~             &                   &                   &                                                   D @@                             ?                                  &                   &                   &                                                     
   @                              @                      @@                              A     è              #DECOMP_INFO               @@                              B                                                                   u #ALLOC_X_REAL C   #ALLOC_X_COMPLEX H   #         @     @X                            C                  #ALLOC_X_REAL%PRESENT D   #VAR E   #OPT_DECOMP F   #OPT_GLOBAL G                 @                            D     PRESENT         D  @                             E                   
               &                   &                   &                                                     
 @@                              F     è             #DECOMP_INFO              
 @@                              G           #         @     @X                            H                  #ALLOC_X_COMPLEX%PRESENT I   #VAR J   #OPT_DECOMP K   #OPT_GLOBAL L                 @                            I     PRESENT         D  @                             J                                  &                   &                   &                                                     
 @@                              K     è             #DECOMP_INFO              
 @@                              L                                                                  u #ALLOC_Y_REAL M   #ALLOC_Y_COMPLEX R   #         @     @X                            M                  #ALLOC_Y_REAL%PRESENT N   #VAR O   #OPT_DECOMP P   #OPT_GLOBAL Q                 @                            N     PRESENT         D  @                             O                   
               &                   &                   &                                                     
 @@                              P     è             #DECOMP_INFO              
 @@                              Q           #         @     @X                            R                  #ALLOC_Y_COMPLEX%PRESENT S   #VAR T   #OPT_DECOMP U   #OPT_GLOBAL V                 @                            S     PRESENT         D  @                             T                                  &                   &                   &                                                     
 @@                              U     è             #DECOMP_INFO              
 @@                              V                                                                  u #ALLOC_Z_REAL W   #ALLOC_Z_COMPLEX \   #         @     @X                            W                  #ALLOC_Z_REAL%PRESENT X   #VAR Y   #OPT_DECOMP Z   #OPT_GLOBAL [                 @                            X     PRESENT         D  @                             Y                   
               &                   &                   &                                                     
 @@                              Z     è             #DECOMP_INFO              
 @@                              [           #         @     @X                            \                  #ALLOC_Z_COMPLEX%PRESENT ]   #VAR ^   #OPT_DECOMP _   #OPT_GLOBAL `                 @                            ]     PRESENT         D  @                             ^                                  &                   &                   &                                                     
 @@                              _     è             #DECOMP_INFO              
 @@                              `                                                        a                                                                                                      b                                                       17                                             c                                                       22           @@                               d                       @@                               e                       @@                               f                       @@                               g                       @@                               h                       @@                               i                       @@                               j                       @@                               k                       @@                               l                       @@                               m                              @               @                'è                   #XST n   #XEN o   #XSZ p   #YST q   #YEN r   #YSZ s   #ZST t   #ZEN u   #ZSZ v   #X1DIST w   #Y1DIST x   #Y2DIST y   #Z2DIST z   #X1CNTS {   #Y1CNTS |   #Y2CNTS }   #Z2CNTS ~   #X1DISP    #Y1DISP    #Y2DISP    #Z2DISP    #X1COUNT    #Y1COUNT    #Y2COUNT    #Z2COUNT    #EVEN                  $                              n                                p          p            p                                        $                              o                               p          p            p                                        $                              p                               p          p            p                                        $                              q            $                   p          p            p                                        $                              r            0                   p          p            p                                        $                              s            <                   p          p            p                                        $                              t            H                   p          p            p                                        $                              u            T                   p          p            p                                        $                              v            `              	     p          p            p                                      $                              w            p              
               &                                                       $                              x            ¸                             &                                                       $                              y                                         &                                                       $                              z            H                            &                                                       $                              {                                        &                                                       $                              |            Ø                            &                                                       $                              }                                         &                                                       $                              ~            h                            &                                                       $                                          °                            &                                                       $                                          ø                            &                                                       $                                          @                            &                                                       $                                                                      &                                                         $                                   Ð                          $                                   Ô                          $                                   Ø                          $                                   Ü                          $                                   à                       @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                                     @                                                         p          p            p                          #         @                                                    #DECOMP_2D_INIT%PRESENT    #NX    #NY    #NZ    #P_ROW    #P_COL    #COMM_LOC    #PERIODIC_BC                  @                                 PRESENT           
  @@                                                   
  @@                                                   
  @@                                                   
   @                                                   
   @                                                   
   @                                                   
 @@                                                  '   p          p            p                          #         @                                                      #         @                                                   #DECOMP_INFO_INIT%ALLOCATED    #DECOMP_INFO_INIT%MAX    #DECOMP_INFO_INIT%MOD    #NX    #NY     #NZ ¡   #DECOMP ¢                 @                                 ALLOCATED               @                                 MAX               @                                 MOD           
  @@                                                   
  @@                                                    
  @@                              ¡                     
D @@                              ¢     è              #DECOMP_INFO    #         @                                 £                   #DECOMP ¤             
D  @                              ¤     è              #DECOMP_INFO    #         @                                 ¥                   #NX ¦   #NY §   #NZ ¨   #PDIM ©   #LSTART ª   #LEND «   #LSIZE ¬             
   @                              ¦                     
   @                              §                     
   @                              ¨                     
   @                              ©                    (   p          p            p                                    D  @                              ª                    )    p          p            p                                    D  @                              «                    *    p          p            p                                    D  @                              ¬                    +    p          p            p                          #         @                                 ­                   #ERRORCODE ®   #MSG ¯             
@ @@                              ®                     
   @                             ¯                    1 #         @                                  °                   #DECOMP ±             D  @                              ±     è              #DECOMP_INFO                  fn#fn    À   ë  b   uapp(DECOMP_2D    «  @   J  MPI %   ë  y       gen@TRANSPOSE_X_TO_Y &   d  ­      TRANSPOSE_X_TO_Y_REAL +     =      TRANSPOSE_X_TO_Y_REAL%SIZE .   N  @      TRANSPOSE_X_TO_Y_REAL%PRESENT *     ¼   a   TRANSPOSE_X_TO_Y_REAL%SRC *   J  ¼   a   TRANSPOSE_X_TO_Y_REAL%DST 1     Y   a   TRANSPOSE_X_TO_Y_REAL%OPT_DECOMP )   _  ³      TRANSPOSE_X_TO_Y_COMPLEX .     =      TRANSPOSE_X_TO_Y_COMPLEX%SIZE 1   O  @      TRANSPOSE_X_TO_Y_COMPLEX%PRESENT -     ¼   a   TRANSPOSE_X_TO_Y_COMPLEX%SRC -   K  ¼   a   TRANSPOSE_X_TO_Y_COMPLEX%DST 4   	  Y   a   TRANSPOSE_X_TO_Y_COMPLEX%OPT_DECOMP %   `	  y       gen@TRANSPOSE_Y_TO_Z &   Ù	  ­      TRANSPOSE_Y_TO_Z_REAL +   
  =      TRANSPOSE_Y_TO_Z_REAL%SIZE .   Ã
  @      TRANSPOSE_Y_TO_Z_REAL%PRESENT *     ¼   a   TRANSPOSE_Y_TO_Z_REAL%SRC *   ¿  ¼   a   TRANSPOSE_Y_TO_Z_REAL%DST 1   {  Y   a   TRANSPOSE_Y_TO_Z_REAL%OPT_DECOMP )   Ô  ³      TRANSPOSE_Y_TO_Z_COMPLEX .     =      TRANSPOSE_Y_TO_Z_COMPLEX%SIZE 1   Ä  @      TRANSPOSE_Y_TO_Z_COMPLEX%PRESENT -     ¼   a   TRANSPOSE_Y_TO_Z_COMPLEX%SRC -   À  ¼   a   TRANSPOSE_Y_TO_Z_COMPLEX%DST 4   |  Y   a   TRANSPOSE_Y_TO_Z_COMPLEX%OPT_DECOMP %   Õ  y       gen@TRANSPOSE_Z_TO_Y &   N  ­      TRANSPOSE_Z_TO_Y_REAL +   û  =      TRANSPOSE_Z_TO_Y_REAL%SIZE .   8  @      TRANSPOSE_Z_TO_Y_REAL%PRESENT *   x  ¼   a   TRANSPOSE_Z_TO_Y_REAL%SRC *   4  ¼   a   TRANSPOSE_Z_TO_Y_REAL%DST 1   ð  Y   a   TRANSPOSE_Z_TO_Y_REAL%OPT_DECOMP )   I  ³      TRANSPOSE_Z_TO_Y_COMPLEX .   ü  =      TRANSPOSE_Z_TO_Y_COMPLEX%SIZE 1   9  @      TRANSPOSE_Z_TO_Y_COMPLEX%PRESENT -   y  ¼   a   TRANSPOSE_Z_TO_Y_COMPLEX%SRC -   5  ¼   a   TRANSPOSE_Z_TO_Y_COMPLEX%DST 4   ñ  Y   a   TRANSPOSE_Z_TO_Y_COMPLEX%OPT_DECOMP %   J  y       gen@TRANSPOSE_Y_TO_X &   Ã  ­      TRANSPOSE_Y_TO_X_REAL +   p  =      TRANSPOSE_Y_TO_X_REAL%SIZE .   ­  @      TRANSPOSE_Y_TO_X_REAL%PRESENT *   í  ¼   a   TRANSPOSE_Y_TO_X_REAL%SRC *   ©  ¼   a   TRANSPOSE_Y_TO_X_REAL%DST 1   e  Y   a   TRANSPOSE_Y_TO_X_REAL%OPT_DECOMP )   ¾  ³      TRANSPOSE_Y_TO_X_COMPLEX .   q  =      TRANSPOSE_Y_TO_X_COMPLEX%SIZE 1   ®  @      TRANSPOSE_Y_TO_X_COMPLEX%PRESENT -   î  ¼   a   TRANSPOSE_Y_TO_X_COMPLEX%SRC -   ª  ¼   a   TRANSPOSE_Y_TO_X_COMPLEX%DST 4   f  Y   a   TRANSPOSE_Y_TO_X_COMPLEX%OPT_DECOMP     ¿  o       gen@UPDATE_HALO !   .  ½      UPDATE_HALO_REAL &   ë  =      UPDATE_HALO_REAL%SIZE )   (  @      UPDATE_HALO_REAL%PRESENT $   h  ¼   a   UPDATE_HALO_REAL%IN %   $  ¼   a   UPDATE_HALO_REAL%OUT '   à  @   a   UPDATE_HALO_REAL%LEVEL ,       Y   a   UPDATE_HALO_REAL%OPT_DECOMP ,   y   @   a   UPDATE_HALO_REAL%OPT_GLOBAL $   ¹   Ã      UPDATE_HALO_COMPLEX )   |!  =      UPDATE_HALO_COMPLEX%SIZE ,   ¹!  @      UPDATE_HALO_COMPLEX%PRESENT '   ù!  ¼   a   UPDATE_HALO_COMPLEX%IN (   µ"  ¼   a   UPDATE_HALO_COMPLEX%OUT *   q#  @   a   UPDATE_HALO_COMPLEX%LEVEL /   ±#  Y   a   UPDATE_HALO_COMPLEX%OPT_DECOMP /   
$  @   a   UPDATE_HALO_COMPLEX%OPT_GLOBAL    J$  g       gen@ALLOC_X    ±$        ALLOC_X_REAL %   <%  @      ALLOC_X_REAL%PRESENT !   |%  ¼   a   ALLOC_X_REAL%VAR (   8&  Y   a   ALLOC_X_REAL%OPT_DECOMP (   &  @   a   ALLOC_X_REAL%OPT_GLOBAL     Ñ&        ALLOC_X_COMPLEX (   _'  @      ALLOC_X_COMPLEX%PRESENT $   '  ¼   a   ALLOC_X_COMPLEX%VAR +   [(  Y   a   ALLOC_X_COMPLEX%OPT_DECOMP +   ´(  @   a   ALLOC_X_COMPLEX%OPT_GLOBAL    ô(  g       gen@ALLOC_Y    [)        ALLOC_Y_REAL %   æ)  @      ALLOC_Y_REAL%PRESENT !   &*  ¼   a   ALLOC_Y_REAL%VAR (   â*  Y   a   ALLOC_Y_REAL%OPT_DECOMP (   ;+  @   a   ALLOC_Y_REAL%OPT_GLOBAL     {+        ALLOC_Y_COMPLEX (   	,  @      ALLOC_Y_COMPLEX%PRESENT $   I,  ¼   a   ALLOC_Y_COMPLEX%VAR +   -  Y   a   ALLOC_Y_COMPLEX%OPT_DECOMP +   ^-  @   a   ALLOC_Y_COMPLEX%OPT_GLOBAL    -  g       gen@ALLOC_Z    .        ALLOC_Z_REAL %   .  @      ALLOC_Z_REAL%PRESENT !   Ð.  ¼   a   ALLOC_Z_REAL%VAR (   /  Y   a   ALLOC_Z_REAL%OPT_DECOMP (   å/  @   a   ALLOC_Z_REAL%OPT_GLOBAL     %0        ALLOC_Z_COMPLEX (   ³0  @      ALLOC_Z_COMPLEX%PRESENT $   ó0  ¼   a   ALLOC_Z_COMPLEX%VAR +   ¯1  Y   a   ALLOC_Z_COMPLEX%OPT_DECOMP +   2  @   a   ALLOC_Z_COMPLEX%OPT_GLOBAL    H2  p       MYTYPE    ¸2  r       REAL_TYPE    *3  r       COMPLEX_TYPE    3  @       MYTYPE_BYTES    Ü3  @       NX_GLOBAL    4  @       NY_GLOBAL    \4  @       NZ_GLOBAL    4  @       COMMLOC    Ü4  @       NRANK    5  @       NPROC &   \5  @       DECOMP_2D_COMM_CART_X &   5  @       DECOMP_2D_COMM_CART_Y &   Ü5  @       DECOMP_2D_COMM_CART_Z    6  o      DECOMP_INFO     7     a   DECOMP_INFO%XST     '8     a   DECOMP_INFO%XEN     Ã8     a   DECOMP_INFO%XSZ     _9     a   DECOMP_INFO%YST     û9     a   DECOMP_INFO%YEN     :     a   DECOMP_INFO%YSZ     3;     a   DECOMP_INFO%ZST     Ï;     a   DECOMP_INFO%ZEN     k<     a   DECOMP_INFO%ZSZ #   =     a   DECOMP_INFO%X1DIST #   =     a   DECOMP_INFO%Y1DIST #   />     a   DECOMP_INFO%Y2DIST #   Ã>     a   DECOMP_INFO%Z2DIST #   W?     a   DECOMP_INFO%X1CNTS #   ë?     a   DECOMP_INFO%Y1CNTS #   @     a   DECOMP_INFO%Y2CNTS #   A     a   DECOMP_INFO%Z2CNTS #   §A     a   DECOMP_INFO%X1DISP #   ;B     a   DECOMP_INFO%Y1DISP #   ÏB     a   DECOMP_INFO%Y2DISP #   cC     a   DECOMP_INFO%Z2DISP $   ÷C  H   a   DECOMP_INFO%X1COUNT $   ?D  H   a   DECOMP_INFO%Y1COUNT $   D  H   a   DECOMP_INFO%Y2COUNT $   ÏD  H   a   DECOMP_INFO%Z2COUNT !   E  H   a   DECOMP_INFO%EVEN    _E         XSTART    óE         XEND    F         XSIZE    G         YSTART    ¯G         YEND    CH         YSIZE    ×H         ZSTART    kI         ZEND    ÿI         ZSIZE    J  ±       DECOMP_2D_INIT '   DK  @      DECOMP_2D_INIT%PRESENT "   K  @   a   DECOMP_2D_INIT%NX "   ÄK  @   a   DECOMP_2D_INIT%NY "   L  @   a   DECOMP_2D_INIT%NZ %   DL  @   a   DECOMP_2D_INIT%P_ROW %   L  @   a   DECOMP_2D_INIT%P_COL (   ÄL  @   a   DECOMP_2D_INIT%COMM_LOC +   M     a   DECOMP_2D_INIT%PERIODIC_BC #   M  H       DECOMP_2D_FINALIZE !   àM  À       DECOMP_INFO_INIT +    N  B      DECOMP_INFO_INIT%ALLOCATED %   âN  <      DECOMP_INFO_INIT%MAX %   O  <      DECOMP_INFO_INIT%MOD $   ZO  @   a   DECOMP_INFO_INIT%NX $   O  @   a   DECOMP_INFO_INIT%NY $   ÚO  @   a   DECOMP_INFO_INIT%NZ (   P  Y   a   DECOMP_INFO_INIT%DECOMP %   sP  T       DECOMP_INFO_FINALIZE ,   ÇP  Y   a   DECOMP_INFO_FINALIZE%DECOMP     Q         PARTITION    «Q  @   a   PARTITION%NX    ëQ  @   a   PARTITION%NY    +R  @   a   PARTITION%NZ    kR     a   PARTITION%PDIM !   ÿR     a   PARTITION%LSTART    S     a   PARTITION%LEND     'T     a   PARTITION%LSIZE     »T  `       DECOMP_2D_ABORT *   U  @   a   DECOMP_2D_ABORT%ERRORCODE $   [U  L   a   DECOMP_2D_ABORT%MSG     §U  T       GET_DECOMP_INFO '   ûU  Y   a   GET_DECOMP_INFO%DECOMP 
#!/bin/bash -l

#FHJ: Logo generated with jp2a utility. See LOGO/donkey2ascii.sh.
#FHJ: Banner was manually generated. It is based on the output of Figlet with font "Varsity"

cat << 'EOF'

                                                                 ..o.          
                                                                .oxxo.         
                                                               .oxxxxo...      
                                                               oxxxxxxo.       
                                                              .oxxxxxxx.       
                                                              .ooooooxxo..     
                                                              .oooooooxo..     
                                                              .oooooxxo...     
                                                       .........oxooo......    
                                                 ............................  
                                             ................................. 
                                          .................................... 
           .          ..oo. ....  .................................oooxxxxxxxo.
     .............oxxxx@ox@@@x@x.....................o...........ooooooooooxx. 
    .o.........oox@x.oo........xxx@@............ooxxxxo..........ooooxxxxxoxo  
    .x........x@xxo...............o@xxo........oxxx@@@xoooooooooooooooxxxo...  
    .o......ox@@o..................oox@o.....ooxxx@xoooxxxxxxxoooooooooooo.... 
    o..ooooo@@xoooo....ooo...........x@o.....ooxxxxo   .oxxxxxxxxxxooooooo.... 
    . .oooo@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....ooooox.     .oxx@@@@@xxoo........ 
      .ooooxxxxxxxooooooxooooooooooooooxo...oooooxx.        ..ox@xxxoo.........
      .ooooooxxxx@xoooooooooooooxxoooooooooooxxxxx.            .oxxooooooooxxo.
     .oooooxxxxx@@@xxxxxxxxxxxxxxxxxxxxxxxxoxxxxo.                .oxxxxxxxxo. 
   ....oooxxxxx@@@xo..oxxx@@@@@@xxxxoxxoooooooxx.                  .oxxxoo..   
  .....ooxxxx@@xo.       ........    .ooooooooxxo                              
  ..oooxxxx@@@o                       .oooooooxoo.                             
  ....oxooxxxxx.                       .ooo..oooo.                             
  .....o.ooxxxxxo.                     .oooooooxo.                             
 ......ooooxxxxxxo.                     .ooooooxoo..                           
........ooxxxxxxxo..                     .o....oxoo...                         
.......ooooxxxxxxxo.                     ........oooo.                         
.ooooooo..ooxxxxoooo.                    .........ooo...                       
..oxo...ooooxxxoooo..                    .ooo......oooo...                     
  .ooooo....o.                            .oxxxoo....ooo....                   
    .oooooo...                              ...ooooo...ooo..                   
       ...                                     .oo.......                      
                                               ....ooo...                      
                      __              __                                       
 ______              [  |            [  |                 ._____  _          _ 
|_   _ \              | |  _          | |                /  ___ \| |        | |
  | |_) | .---.  _. _.| | / |  .---.  | | .---.  _    _ / /    \_|\ \  /\  / / 
  |  __'./ /__\\[ /`\_| '' <  / /__\\ | |/ /__\\| \  | | |   _____ \ \/  \/ /  
 _| |__| | \__. | |   | |`\ \ | \___. | || \___. \ \/ / \ \.___| |  \  /\  /   
|_______/ \.__./[_]  [__|  \_] \.__./[___]\.__./  \  /   \.____./    \/  \/    
                                                  / /                          
                                                 /_/                           

EOF

if [ "$1" = ".real" ]; then
	echo -e "\n*************************     Building REAL flavor    **************************\n"
elif [ "$1" = ".cplx" ]; then
	echo -e "\n************************    Building COMPLEX flavor     ************************\n"
fi

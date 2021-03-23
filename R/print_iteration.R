print_iteration <-
function(para, p, q, diff, dall, method, finish, t, m){
            
 if(method=="Delta_l_x"){ 
             if(p==1){ if(q==1){ cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'Delta_lx=',
                                   diff, 'Stop?=', finish,  "m=", m,'\n') }
                        if(q==2){ cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Delta_lx=',
                                        diff, 'Stop?=', finish,  "m=", m,'\n')  }
                        }
            if(p==2){ if(q==1){ cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'Delta_lx=',
                                       diff, 'Stop?=', finish,  "m=", m,'\n')  }
                        if(q==2){ cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Delta_lx=',
                                        diff, 'Stop?=', finish,  "m=", m,'\n')  }
                        }
            
 }



if(method=="MA_1"){           
            if(p==1){ if(q==1){cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'Max parameter diff=',
                                        diff, 'Stop?=', finish,  "m=", m,'\n')  }
                        if(q==2){cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Max parameter diff=',
                                        diff, 'Stop?=', finish,  "m=", m,'\n')  }
            }
            
            if(p==2){ if(q==1){cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'Max parameter diff=',
                                  diff, 'Stop?=', finish,  "m=", m,'\n')  }
                      if(q==2){ cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Max parameter diff=',
                                   diff, 'Stop?=', finish,  "m=", m,'\n') }
            }
}
if(method=="MA_k"){           
            if(p==1){   if(q==1){cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'Max parameter diff=',
                                        diff, 'Stop?=', finish,  "m=", m,'\n')}                       
                        if(q==2){cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Max parameter diff=',
                                        diff, 'Stop?=', finish,  "m=", m,'\n')}
                     }
            
            if(p==2){   if(q==1){cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'Max parameter diff=',
                                 diff, 'Stop?=', finish,  "m=", m,'\n') }
                        if(q==2){cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'Max parameter diff=',
                                 diff, 'Stop?=', finish,  "m=", m,'\n') }
            }
}



if(method=="ALL"){           
            if(p==1){   if(q==1){cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'iteration diff=',
                                     dall, 'Stop?=', finish,  "m=", m,'\n')}                       
                        if(q==2){cat(t-1, 'phi1=', para[1],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'iteration diff=',
                                     dall, 'Stop?=', finish,  "m=", m,'\n')}
            }
            
            if(p==2){   if(q==1){cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'iteration diff=',
                                     dall, 'Stop?=', finish,  "m=", m,'\n') }
                        if(q==2){cat(t-1, 'phi1=', para[1], 'phi2=', para[2],  'theta1=', para[(p+1)], 'theta2=', para[(p+2)], 'iteration diff=',
                                     dall, 'Stop?=', finish,  "m=", m,'\n') }
            }
}


}

# functions for converting sensor permittivity [-] to voltage  [Volts] and back

# voltage to permittivity ####
V2eps = function(V, type, temp=NULL)
{  
  
  #PR2
  #default equation according to manual (PR2_user_manual_version_5.0.pdf, eq. 2)
  # (PR2_SDI-12-_User_Manual_version_4_1.pdf, eq. 2)
  if (type=="PR2")
    return(c(eps =(1.125 - 5.53*V + 67.17*V^2 - 234.42*V^3 + 413.56*V^4 - 356.68*V^5 + 121.53*V^6)^2)) 
  
  #Theta-probe, polynomial
  #convert to epsilon, eq. 1 of Theta Probe user manual, p.12
  if (type=="Theta Probe polynomial")
    return(c(eps = (1.07 + 6.4*V-6.4*V^2+4.7*V^3 )^2))
  
  #theta-probe, lookup values reconstructed from table p. 14
  #converted to epsilon, eq. 1 of Theta Probe user manual, p.12
  #this is more general than "Theta Probe polynomial", as it also fits above 0.55 V
  if (type=="Theta Probe")
    return(c(eps = .V2eps_Theta_Probe_table(V)))
    #return(c(eps = do.call(globvars$V2eps_Theta_Probe_table, args=list(V), envir = globvars))) #geht!
  
    #return(c(eps = V))
  
  if(type=="picoSMS_linear" || type=="picoSMS"){
    if (missing(temp) || is.null(temp)) {
      stop("'temp' must be provided and cannot be NULL.")
    }
    if(type=="picoSMS"){
      warning("Using linear function for picoSMS. State type == 'picoSMS_poly' for polynomial function")
    }
    return(c(eps=2.3811 + 70.0811 * V - 0.0337 * temp))
  }
  
  if(type=="picoSMS_poly"){
    if (missing(temp) || is.null(temp)) {
      stop("'temp' must be provided and cannot be NULL.")
    }
    return(c(eps=1.9284 + 43.9328 * V + 0.0181 * temp + 26.4293 * V^2 - 0.0942 * V * temp - 0.0005 * temp^2))
  }
  
  stop("type must be 'PR2', 'Theta Probe', 'Theta Probe polynomial', 'picoSMS_linear' or 'picoSMS_poly'.")
}



#permittivity to voltage ####
  .eps2V = list() #internal list containing functions for each sensor type
  #find inverse functions: since there are no easy analytical inverse, we use piecewise linear inverse functions 
  for (probe_type in c("PR2", "Theta Probe polynomial")) #for "Theta Probe", this is done later in .onLoad()
  {  
    minV = optimize(V2eps, interval=c(-0.2, 0.1), type = probe_type)$minimum  #find minimum of regression, i.e. start of monotonically increasing parabola
  
    V_synth = seq(from=minV, to=1.35, length.out=500) #range of Voltage values [V] used for lookup table
    eps_synth = V2eps(V_synth, type=probe_type) #range of epsilon values
    .eps2V[[probe_type]] = approxfun(x=eps_synth, y=V_synth) #use piecewise linear approximation

    # # illustration ####
    # V2 =.eps2V[[probe_type]](eps_synth) #approximate inverse of default regression
    # plot(V_synth, eps_synth, type="l", lwd=2, main=probe_type) #, xlim=c(0,0.1), ylim=c(0,))
    # lines(V2, eps_synth, col="red", lty="dashed")
  }
  
  eps2V = function(eps, type, temp = NULL) #wrapper for internal function created above
  {  
    
    if(type == "picoSMS_linear" || type == "picoSMS_poly"|| type == "picoSMS"){
      
      if(type == "picoSMS_linear"|| type == "picoSMS"){
        if(type=="picoSMS"){
          warning("Using linear function for picoSMS. State type == 'picoSMS_poly' for polynomial function")}
        return(V = (eps - 2.3811 + 0.0337 * temp)/70.0811)
      }
      
      if(type == "picoSMS_poly"){
        a= 26.4293
        b= 43.9328 - 0.0942 * temp
        c= 1.9284 + 0.0181 * temp - 0.0005 * temp^2 - eps
        root_term = b^2 - 4*a*c
        return((-b + sqrt(root_term))/(2*a)) # negative is not valid since it leads to negative V
      }
    }else{
    
    #replace e.g. "PR2, analogue" by "PR2"
      type = sub(x = type, pattern = paste0("^(", paste0(names(.eps2V), collapse="|" ), ").*"), replacement ="\\1")
      if (!(type %in% names(.eps2V)))
        stop("Unknown probe type. Must be one of ('", paste0(names(.eps2V), collapse="', '"),"').")
      return(.eps2V[[type]](eps))
    }
    
  }


  
module hcore
    use variables 
    use physical_properties 

    contains
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function h_coefficient(t,n_axis)                 !!!!t代表时间t，n_axis代表轴向划分网格数
    integer t
    if (module_identification==1) then
        if (coolant_kind==1) then
            Re=(11096.0-1.3236*Temperature(t,n_axis,n_radial+3))*coolant_speed*De/(4.94e-4*exp(757.1/&
                &Temperature(t,n_axis,n_radial+3)))
            Pr=(159.0-0.0272*Temperature(t,n_axis,n_radial+3)+7.12e-6*Temperature(t,n_axis,n_radial+3&
                &)**2)*(4.94e-4*exp(757.1/Temperature(t,n_axis,n_radial+3)))/(3.61+0.01517*Temperatur&
                &e(t,n_axis,n_radial+3))
            Pe=Re*Pr
            Nu=7.55*(pitch/d(t,n_axis,n_radial+2))-20.*((pitch/d(t,n_axis,n_radial+2))**-13)+3.67/90.&
                &*((pitch/d(t,n_axis,n_radial+2))**-2)*(Pe**(0.56+0.19*(pitch/d(t,n_axis,n_radial+2))))
            h_coefficient=Nu*(3.61+0.01517*Temperature(t,n_axis,n_radial+3))/De
        end if
    else
            Re=(11096.0-1.3236*Temperature_transient(t,n_axis,n_radial+3))*coolant_speed*De/(4.94e-4*exp(757.1/&
                &Temperature_transient(t,n_axis,n_radial+3)))
            Pr=(159.0-0.0272*Temperature_transient(t,n_axis,n_radial+3)+7.12e-6*Temperature_transient(t,n_axis,n_radial+3&
                &)**2)*(4.94e-4*exp(757.1/Temperature_transient(t,n_axis,n_radial+3)))/(3.61+0.01517*Temperatur&
                &e_transient(t,n_axis,n_radial+3))
            Pe=Re*Pr
            Nu=7.55*(pitch/d(t,n_axis,n_radial+2))-20.*((pitch/d(t,n_axis,n_radial+2))**-13)+3.67/90.&
                &*((pitch/d(t,n_axis,n_radial+2))**-2)*(Pe**(0.56+0.19*(pitch/d(t,n_axis,n_radial+2))))
            h_coefficient=Nu*(3.61+0.01517*Temperature_transient(t,n_axis,n_radial+3))/De
    end if
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!对流换热系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Kc(i,j)   
    real(8)::i,j,Kc
     Kc=(27.67 + 1.39d-2*((i+j)/2.) - 2.42d-5*((i+j)/2.)**2 + 8.29d-9*((i+j)/2.)**3)
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!包壳热导率!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function Kg(coolant_kind)
    integer coolant_kind
        if(coolant_kind==1) then
            Kg=5678
        end if
    end function
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!气隙热导率Kg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Bu_calculation(time1,n_a)
    integer n_a,time1,n_r,day
    real(8)::sum_q
    day=1
    sum_q=0.
    Bu_average=4.*p_line(n_a)/(pi*(d(time1,n_a,n_radial)**2))*day/(UO2_density*238./270.)/10.**6.
    Bu(time1,n_a,:)=Bu_average
    
    if (time1==1) then
        Bu(time1,n_a,1)=Bu(time1,n_a,n_radial-2)*(1.-Bu(time1,n_a,n_radial-2)/300.)
        sum_q=sum_q+Bu(time1,n_a,1)*pi*(0.1**2.)
        do n_r=2,n_radial-2
            Bu(time1,n_a,n_r)=Bu(time1,n_a,n_r-1)+(Bu(time1,n_a,n_radial-2)-Bu(time1,n_a,1))/8.
            sum_q=sum_q+Bu(time1,n_a,n_r)*pi*((0.1*n_r)**2.-(0.1*(n_r-1))**2.)
        end do
        Bu(time1,n_a,n_radial-1)=(Bu(time1,n_a,n_radial-2)*pi*(1.0**2.)-sum_q)/(pi*(1.-0.9**2.))
        Bu(time1,n_a,n_radial)=Bu(time1,n_a,n_radial-1)
    else
        do i=1,n_radial
            Bu_old(time1,n_a,i)=Bu(time1-1,n_a,i)
        end do
        Bu(time1,n_a,1)=Bu(time1,n_a,n_radial-2)*(1.-Bu(time1,n_a,n_radial-2)/300.)
        sum_q=sum_q+Bu(time1,n_a,1)*pi*(0.1**2.)
        do n_r=2,n_radial-2
            Bu(time1,n_a,n_r)=Bu(time1,n_a,n_r-1)+(Bu(time1,n_a,n_radial-2)-Bu(time1,n_a,1))/8.
            sum_q=sum_q+Bu(time1,n_a,n_r)*pi*((0.1*n_r)**2.-(0.1*(n_r-1))**2.)
        end do
        Bu(time1,n_a,n_radial-1)=(Bu(time1,n_a,n_radial-2)*pi*(1.0**2.)-sum_q)/(pi*(1.-0.9**2.))
        Bu(time1,n_a,n_radial)=Bu(time1,n_a,n_radial-1)
        do i=1,n_radial
            Bu(time1,n_a,i)=Bu(time1,n_a,i)+Bu_old(time1,n_a,i)
        end do
        
    end if
    end subroutine
    
    
    subroutine P_S_calculation(time1,n_a)
    integer time1,n_a
    real(8)::P_SQUARE_AVERAGE,sum_p
    sum_p=0.
    p_square_average=4.*P_line(n_a)/(pi*(d(time1,n_a,n_radial)**2))
    P_SQUARE(:,n_a,:)=p_square_average
    P_SQUARE(time1,n_a,1)=P_SQUARE(time1,n_a,n_radial-2)*(1-Bu_average/300.)
    sum_p=sum_p+P_SQUARE(time1,n_a,1)*pi*(0.1**2.)
    do n_r=2,n_radial-2
        P_SQUARE(time1,n_a,n_r)=P_SQUARE(time1,n_a,n_r-1)+(P_SQUARE(time1,n_a,n_radial-2)-P_SQUARE(time1,n_a,1))/8.
        sum_p=sum_p+P_SQUARE(time1,n_a,n_r)*pi*((0.1*n_r)**2.-(0.1*(n_r-1))**2.)
    end do
    P_SQUARE(time1,n_a,n_radial-1)=(P_SQUARE(time1,n_a,n_radial-2)*pi*(1.0**2.)-sum_p)/(pi*(1.-0.9**2.))
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!燃耗与体积功率计算Bu&P!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fuel_coefficient_calculation(time1,n_a,n_r)
    integer time1,n_a,n_r
    real(8)::S,P,B,k0,fd,fp,fm,fr
    S=1.5           !!!!!!!!!!形状因子
    P=0.05          !!!!!!!!!!芯块孔隙率
    if (module_identification==1) then
        k0=1./(0.0375+(2.165d-4)*Temperature(time1,n_a,n_r))+4.715d9*exp(-16361./Temperature(time1,n_a,n_r))/Temperature(time1,n_a,n_r)**2
        !!!!!!!!!!!!!!!!k0-未考虑燃耗作用的二氧化铀导热系数
        FD=(1.09/(Bu(time1,n_a,n_r)/9.383)**3.265+0.0643/((Bu(time1,n_a,n_r)/9.383)**0.5)*(Temperature(time1,n_a,n_r)**0.5)&
         &)*ATAN2(1.0,(1.09/((Bu(time1,n_a,n_r)/9.383)**3.265)+0.0643/((Bu(time1,n_a,n_r)/9.383)**0.5)*(Temperature(time1,n_a,n_r)**0.5)))
        !!!!!!!!!!!!!!!!FD-溶解因子
        FP=1+(0.019*(Bu(time1,n_a,n_r)/9.383)/(3.-0.019*(Bu(time1,n_a,n_r)/9.383)))*(1./(1.+EXP(-(Temperature(time1,n_a,n_r)-1200)/100.)))
        !!!!!!!!!!!!!!!!FP-沉淀因子
        FM=(1-P)/(1+(S-1)*P)!!!!!!!!!FM-孔隙因子
        FR=1.-0.2/(1.+EXP((Temperature(time1,n_a,n_r)-900.)/80.))!!!!!!!!!!!!!FR辐照因子
        k_fuel(time1,n_a,n_r)=k0*FD*FM*FR               !!!!!!!!二氧化铀真实导热系数
    else
        k0=1./(0.0375+(2.165d-4)*Temperature_transient(time1,n_a,n_r))+4.715d9*exp(-16361./Temperature_transient(time1,n_a,n_r))/Temperature_transient(time1,n_a,n_r)**2
        !!!!!!!!!!!!!!!!k0-未考虑燃耗作用的二氧化铀导热系数
        FD=(1.09/(Bu(time,n_a,n_r)/9.383)**3.265+0.0643/((Bu(time,n_a,n_r)/9.383)**0.5)*(Temperature_transient(time1,n_a,n_r)**0.5)&
         &)*ATAN2(1.0,(1.09/((Bu(time,n_a,n_r)/9.383)**3.265)+0.0643/((Bu(time,n_a,n_r)/9.383)**0.5)*(Temperature_transient(time1,n_a,n_r)**0.5)))
        !!!!!!!!!!!!!!!!FD-溶解因子
        FP=1+(0.019*(Bu(time,n_a,n_r)/9.383)/(3.-0.019*(Bu(time,n_a,n_r)/9.383)))*(1./(1.+EXP(-(Temperature_transient(time1,n_a,n_r)-1200)/100.)))
        !!!!!!!!!!!!!!!!FP-沉淀因子
        FM=(1-P)/(1+(S-1)*P)!!!!!!!!!FM-孔隙因子
        FR=1.-0.2/(1.+EXP((Temperature_transient(time1,n_a,n_r)-900.)/80.))!!!!!!!!!!!!!FR辐照因子
        k_fuel(time1,n_a,n_r)=k0*FD*FM*FR               !!!!!!!!二氧化铀真实导热系数   
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!UO2芯块热导率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
end module
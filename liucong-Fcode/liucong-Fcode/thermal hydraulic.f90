    module thermal_hydraulic
    use variables
    use physical_properties
    use hcore
    use match
    contains



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine power_calculation
    real(8)::a,b

    p_line_max=p_line_average*p_line_factor
    x=0.                                !!!!假定功率函数截尾位置
    Fq=0.                               !!!!迭代赋值
    do while(abs(Fq-p_line_factor)>0.005)
        sin_average=(cos(x)-cos(pi-x))/(pi-2*x)
        Fq=1./sin_average
        if (Fq-p_line_factor>0) then
            x=x+0.001
        else
            x=x-0.001
        end if
    end do
    do i=1,n_axis
        p_line(i)=p_line_max*sin(x+(i-1)*(pi-2*x)/(1.0*(n_axis-1)))
    end do
    do i=1,n_axis
    a=0
    b=0
    a=(i-1)*cladding_length/((n_axis-1)*1.0)
    b=(i-2)*cladding_length/((n_axis-2)*1.0)
    Q_line(i)=p_line_max*cladding_length/(pi-2*x)*(cos(x+b*(pi-2*x)/cladding_length)-cos(x+a*(pi-2*x)/cladding_length))
    end do
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!功率计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine coolant_T_calculation(n_a)
    real(8)::a,b
    if (module_identification==1)then
        a=0
        a=(n_a-1)*cladding_length/((n_axis-1)*1.0)
        Q_line(n_a)=p_line_max*cladding_length/(pi-2*x)*(cos(x)-cos(x+a*(pi-2*x)/cladding_length))
        coolant_T=Q_line(n_a)/coolant_M_flow/coolant_Cp+coolant_T_in
    else
        if (n_a==1) then
            Q_line(n_a)=0.
            coolant_T_transient(n_a)=coolant_T_in
        else
            coolant_T_transient(n_a)=coolant_T_transient(n_a-1)+Q_line(n_a)/coolant_M_flow/coolant_Cp
        end if
    end if
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!冷却剂温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine axis_T_calculation(time,n_a)
    integer time
    cladding_T_surface=Temperature(time,n_a,n_radial+3)+p_line(n_a)&
        &/(pi*d(time,n_a,n_radial+2)*h_coefficient(time,n_a))
    cladding_T=cladding_T_surface                                       !!!!!!假设包壳内表面温度
    cladding_T_internal=cladding_T_surface+P_line(n_a)*log(&            !!!!!包壳内表面初始迭代值
    &d(time,n_a,n_radial+2)/(d(time,n_a,n_radial+2)-2*cl&
        &adding_width))/(2*pi*Kc(cladding_T,cladding_T_surface))
    Tc1=cladding_T_internal-cladding_T                                  !!!!!初始温差
    Tc=0                                                                !!!!!!假设温差
    !包壳迭代计算开始
    do while((Tc1-Tc)/Tc1>0.005)
        Tc=Tc1
        cladding_T_internal=cladding_T_surface+p_line(n_a)&
            &*log(d(time,n_a,n_radial+2)/(d(time,n_a,n_radial+2)&
            &-2*cladding_width))/(2*pi*Kc(cladding_T,cladding_T_surface))
        Tc1=abs(cladding_T-cladding_T_internal)
        cladding_T=cladding_T_internal
    enddo
    !包壳迭代计算结束
    pellet_T_surface=cladding_T_internal+P_line_max*log((d(time,n_a&!气隙换热，计算芯块外表面温度，Kg为假定
    &,n_radial+2)-2*cladding_width)/pellet_diame&
        &ter)/(2*pi*Kg(coolant_kind))
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!轴向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度更新!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度更新!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度更新!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine axis_T_update(time,n_a)
    integer time
    Temperature(time,n_a,n_radial+2)=cladding_T_surface
    Temperature(time,n_a,n_radial+1)=cladding_T_internal
    do i=1,n_radial
        Temperature(time,n_a,i)=pellet_T_surface
    end do
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度更新!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度更新!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度更新!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!径向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!径向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!径向温度计算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fuel_t_calculation(time1,n_a)
    integer time1,n_a
    dimension A(n_radial,n_radial),B(n_radial),X(n_radial),Tc(n_radial)
    double precision A,B,X
    real(8)::tt
    X(:)=0                                                  !!!!矩阵初始化
    tt=100.
    A(:,:)=0.
    do while(tt>0.01)
        do i=1,n_radial
            call fuel_coefficient_calculation(time1,n_a,i)
        end do
        !!!!!!!系数矩阵赋值开始!!!!!!!!
        A(1,1)=1.
        A(1,2)=-1.
        B(1)=0.5*P_SQUARE(TIME1,n_a,1)*(d(time1,n_a,2)**2)/4./k_fuel(time1,n_a,1)
        do i=2,n_radial-1
            A(i,i-1)=(d(time1,n_a,i)+d(time1,n_a,i-1))/(d(time1,n_a,i)-d(time1,n_a,i-1))/2.*k_fuel(time1,n_a,i-1)
            A(i,i+1)=(d(time1,n_a,i)+d(time1,n_a,i+1))/(d(time1,n_a,i+1)-d(time1,n_a,i))/2.*k_fuel(time1,n_a,i)
            A(i,i)=-(A(i,i-1)+A(i,i+1))
            B(i)=-1.*P_SQUARE(TIME1,n_a,i)*d(time1,n_a,i)/2.*(d(time1,n_a,i+1)-d(time1,n_a,i-1))/4.
        end do
        A(n_radial,n_radial)=1.
        B(n_radial)=pellet_T_surface
        !!!!!!!系数矩阵赋值结束!!!!!!!!
        call TDMA(n_radial,A,B,X)                            !!!!!TDMA计算，后续继续采取该计算迭代芯块温度值
        do i=1,n_radial
            Tc(i)=x(i)-Temperature(time1,n_a,i)
            Temperature(time1,n_a,i)=x(i)
        end do
        tt=maxval(abs(Tc))
    end do
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!径向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!径向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!径向温度计算结束!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度计算瞬态!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度计算瞬态!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度计算瞬态!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine T_transient_calculation(time11,n_a)
    integer i,time11,n_a
    dimension A(n_radial+2,n_radial+2),B(n_radial+2),X(n_radial+2),Tc(n_radial+2),transient(n_radial+2)
    double precision A,B,X,transient
    real(8)::tt,rw,re,d_rw,d_re,heat
    cladding_T_surface=Temperature_transient(time11,n_a,n_radial+3)+p_line(n_a)&
        &/(pi*d(time11,n_a,n_radial+2)*h_coefficient(time11,n_a))
    do i=1,n_radial+2
        Temperature_transient(time11,n_a,i)=cladding_T_surface
    end do








    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    X(:)=0                                                  !!!!矩阵初始化
    tt=100.
    A(:,:)=0.
    do while(tt>0.01)
        do i=1,n_radial
            call fuel_coefficient_calculation(time11,n_a,i)
        end do
        !!!!!!!系数矩阵赋值开始!!!!!!!!
        if (time11==1) then
            transient(1)=UO2_denstiy*UO2_Cp*(d(time11,n_a,2)/4.)**2./AT/2.
            A(1,1)=-(1.*k_fuel(time11,n_a,1)/2.+transient(1))
            A(1,2)=1.*k_fuel(time11,n_a,1)/2.
            B(1)=-(0.5*P_SQUARE(time11,n_a,1)*(d(time11,n_a,2)/4.)**2.+transient(1)*Temperature(time,n_a,1))
            do i=2,n_radial+1
                rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
                re=(d(time11,n_a,i)+d(time11,n_a,i+1))/4.
                d_rw=(d(time11,n_a,i)-d(time11,n_a,i-1))/2.
                d_re=(d(time11,n_a,i+1)-d(time11,n_a,i))/2.
                if (i<n_radial) then
                    transient(i)=UO2_density*UO2_Cp/AT*(re-rw)*(re+rw)/2.
                    A(i,i-1)=rw/d_rw*k_fuel(time11,n_a,i-1)
                    A(i,i+1)=re/d_re*k_fuel(time11,n_a,i)
                    A(i,i)=-(A(i,i-1)+A(i,i+1)+transient(i))
                    B(i)=-1.*(P_SQUARE(time11,n_a,i)*(re-rw)*(re+rw)/2.+transient(i)*Temperature(time,n_a,i))
                else if (i==n_radial) then
                    A(i,i-1)=k_fuel(time11,n_a,i)/d_re
                    A(i,i+1)=Kg(coolant_kind)
                    A(i,i)=-(A(i,i-1)+A(i,i+1))
                    B(i)=0
                else if (i==n_radial+1) then
                    A(i,i-1)=Kg(coolant_kind)*d(time11,n_a,n_radial)/d(time11,n_a,n_radial+1)
                    A(i,i+1)=Kc(Temperature_transient(time11,n_a,i),cladding_T_surface)/d_rw
                    A(i,i)=-(A(i,i-1)+A(i,i+1))
                    B(i)=0
                end if

            end do
            A(n_radial+2,n_radial+2)=1.
            B(n_radial+2)=cladding_T_surface
        else
            transient(1)=UO2_denstiy*UO2_Cp*(d(time11,n_a,2)/4.)**2./AT/2.
            A(1,1)=-(1.*k_fuel(time11,n_a,1)/2.+transient(1))
            A(1,2)=1.*k_fuel(time11,n_a,1)/2.
            B(1)=-(0.5*P_SQUARE(time11,n_a,1)*(d(time11,n_a,2)/4.)**2.+transient(1)*Temperature_transient(time11-1,n_a,1))
            do i=2,n_radial+1
                rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
                re=(d(time11,n_a,i)+d(time11,n_a,i+1))/4.
                d_rw=(d(time11,n_a,i)-d(time11,n_a,i-1))/2.
                d_re=(d(time11,n_a,i+1)-d(time11,n_a,i))/2.
                if (i<n_radial) then
                    transient(i)=UO2_density*UO2_Cp/AT*(re-rw)*(re+rw)/2.
                    A(i,i-1)=rw/d_rw*k_fuel(time11,n_a,i-1)
                    A(i,i+1)=re/d_re*k_fuel(time11,n_a,i)
                    A(i,i)=-(A(i,i-1)+A(i,i+1)+transient(i))
                    B(i)=-1.*(P_SQUARE(time11,n_a,i)*(re-rw)*(re+rw)/2.+transient(i)*Temperature_transient(time11-1,n_a,i))
                else if (i==n_radial) then
                    A(i,i-1)=k_fuel(time11,n_a,i)/d_re
                    A(i,i+1)=Kg(coolant_kind)
                    A(i,i)=-(A(i,i-1)+A(i,i+1))
                    B(i)=0
                else if (i==n_radial+1) then
                    A(i,i-1)=Kg(coolant_kind)*d(time11,n_a,n_radial)/d(time11,n_a,n_radial+1)
                    A(i,i+1)=Kc(Temperature_transient(time11,n_a,i),cladding_T_surface)/d_rw
                    A(i,i)=-(A(i,i-1)+A(i,i+1))
                    B(i)=0
                end if

            end do
            A(n_radial+2,n_radial+2)=1.
            B(n_radial+2)=cladding_T_surface
        end if
        !!!!!!!系数矩阵赋值结束!!!!!!!!
        call TDMA(n_radial+2,A,B,X)                            !!!!!TDMA计算，后续继续采取该计算迭代芯块温度值
        do i=1,n_radial+2
            Tc(i)=x(i)-Temperature_transient(time11,n_a,i)
            Temperature_transient(time11,n_a,i)=x(i)
        end do
        tt=maxval(abs(Tc))
    end do







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!计算储热
    Tc(:)=0.
    heat=0.
    if (time11==1) then
        do i=1,n_radial+2
            Tc(i)=Temperature_transient(time11,n_a,i)-Temperature(time,n_a,i)
            if (i==1) then
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*((d(time11,n_a,i)/2.)**2)/2.
            else if(i<n_radial) then
                rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
                re=(d(time11,n_a,i)+d(time11,n_a,i+1))/4.
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(re**2-rw**2)
            else if(i==n_radial) then
                rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(d(time11,n_a,i)**2/4.-rw**2)
            else
                rw=d(time11,n_a,n_radial+1)/2.
                re=d(time11,n_a,n_radial+2)/2.
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(re**2-rw**2)/2.
            end if
        end do
    else
        do i=1,n_radial+2
            Tc(i)=Temperature_transient(time11,n_a,i)-Temperature_transient(time11-1,n_a,i)
            if (i==1) then
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*((d(time11,n_a,i)/2.)**2)/2.
            else if(i<n_radial) then
                rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
                re=(d(time11,n_a,i)+d(time11,n_a,i+1))/4.
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(re**2-rw**2)
            else if(i==n_radial) then
                rw=(d(time11,n_a,i)+d(time11,n_a,i-1))/4.
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(d(time11,n_a,i)**2/4.-rw**2)
            else
                rw=d(time11,n_a,n_radial+1)/2.
                re=d(time11,n_a,n_radial+2)/2.
                heat=heat+Tc(i)*UO2_Cp*UO2_density*pi*(re**2-rw**2)/2.
            end if
        end do
    end if


    do i=2,n_axis
    Q_line(i)=Q_line(i)-heat*length_every        !热量更新
    end do
    heat_dv=heat-heat_old                       !热量差值

    heat_old=heat




    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度计算瞬态!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度计算瞬态!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!温度计算瞬态!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module
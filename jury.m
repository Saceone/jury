function jury(v)

%jury(v)
%   v  = symbolic characteristic polynomial of the system to be studied

n=length(v);
digits(5);

if strcmp(class(v),'double')
    inestable=0;
    fprintf('\nAplicamos el criterio de Jury al polinomio %s\n',char(vpa(poly2sym(v,'z'))))
    if n>=3
        M=zeros((2*n-4),n);
        M(2,:)=v(n:-1:1);
        M(1,:)=v(1:1:n);
        
        for i=1:n-3
            for j=1:n-i
               M(2*i+1,j)=-det( [M(2*i-1,1) M(2*i-1,j+1);M(2*i,1) M(2*i,j+1)] );
               M(2*i+2,n-i+1-j)=M(2*i+1,j); 
            end
        end

        fprintf('\n')
        display('Construimos la tabla de Jury:')
        fprintf('\n')
        TABLA_JURY=vpa(M);
        pretty(TABLA_JURY)
    display('Comprobamos las condiciones del criterio:')
    izq=char(vpa(M(2,1)));
    der=char(vpa(M(1,1)));
    if abs(M(2,1))<abs(M(1,1))
        val='OK';
    else
        val='NO';
        inestable=1;
    end
    fprintf('------> |%s| < |%s| ==> %s\n',izq,der,val);
    i=3;
    while i<length(M)
        izq=char(vpa(M(i+1,1)));
        der=char(vpa(M(i,1)));
        if abs(M(i+1,1))>abs(M(i,1))
            val='OK';
        else
            val='NO';
            inestable=1;
        end
        fprintf('------> |%s| > |%s| ==> %s\n',izq,der,val);
        i=i+2;
    end
    else
        display('El polinomio característico es de orden menor que 3: no ha lugar construir la tabla')
        display('Comprobamos las condiciones del criterio:')
        izq=char(vpa((1)));
        der=char(vpa(v(n)));
        if abs(v(n))<abs(v(1))
            val='OK';
        else
            val='NO';
            inestable=1;
        end
        fprintf('------> |%s| > |%s| ==> %s\n',izq,der,val);
    end
    if inestable==1
        fprintf('\nSISTEMA INESTABLE\n\n')
    else
        fprintf('\nSISTEMA ESTABLE\n\n')
    end
    display('¡Gracias por usar la función jury!')
else
    fprintf('\nAplicamos el criterio de Jury al polinomio %s\n\n',char(vpa(poly2sym(v,'z'))))
    %%%%%%%%%%CONDICIONES PREVIAS%%%%%%%%%%%%%%%%%

    %%%P(1)>0%%%
    syms z k
    eq=poly2sym(v,'z');
    eq=subs(eq,z,1);
    display('De la condición P(1)>0 se deduce que:')
    intervalo_1=vpa(solve(eq>0,k))
    %%%%%%%%%%%%

    %%%P(-1)<>0 según orden%%%
    if mod(n-1,2) == 0 %orden par
        eq1=poly2sym(v,'z');
        eq1=subs(eq1,z,-1);
        display('De la condición P(-1)>0 (por ser p(z) de orden par) se deduce que:')
        intervalo_2=vpa(solve(eq1>0,k))
    else %orden impar
        eq1=poly2sym(v,'z');
        eq1=subs(eq1,z,-1);
        display('De la condición P(-1)<0 (por ser p(z) de orden impar) se deduce que:')
        intervalo_2=vpa(solve(eq1<0,k))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    w = warning ('off','all');
    display('Por tanto, un primer intervalo es la intersección de los anteriores: ')
    if mod(n-1,2) == 0
        intervalo_previo=vpa(solve(eq>0,eq1>0))
    else
        intervalo_previo=vpa(solve(eq>0,eq1<0))
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%CONSTRUIMOS LA TABLA%%%%%%%%%%%%%%%
    if n>3
        M=sym(zeros((2*n-4),n));
        M(2,:)=v(n:-1:1);
        M(1,:)=v(1:1:n);
        
        for i=1:n-3
            for j=1:n-i
               M(2*i+1,j)=-det( [M(2*i-1,1) M(2*i-1,j+1);M(2*i,1) M(2*i,j+1)] );
               M(2*i+2,n-i+1-j)=M(2*i+1,j); 
            end
        end
        
        fprintf('\n')
        display('Construimos la tabla de Jury:')
        fprintf('\n')
        TABLA_JURY=vpa(M);
        pretty(TABLA_JURY)
        %%%%%%%%%%%COMPROBACIONES FINALES%%%%%%%%%%%%%
        display('Para determinar por completo el intervalo, imponer:')
        izq=char(vpa(M(2,1)));
        der=char(vpa(M(1,1)));
        fprintf('------> |%s| < |%s|\n',izq,der);
        i=3;
        while i<length(M)
            izq=char(vpa(M(i+1,1)));
            der=char(vpa(M(i,1)));
            fprintf('------> |%s| > |%s|\n',izq,der);
            i=i+2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        display('El polinomio característico es de orden menor que 3: no ha lugar construir la tabla')
        display('Para determinar por completo el intervalo, imponer:')
        izq=char(vpa(v(n)));
        der=char(vpa(v(1)));
        fprintf('------> |%s| < |%s|\n',izq,der);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n\n');
    display('¡Gracias por usar la función jury!')
end
end

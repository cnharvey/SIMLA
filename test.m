clear

clc

gofactor=30;
stopfactor=20;
skipfactor=0.1;

counter=0;



while(true)
    counter=counter+1;
       
    
    
    if counter==5
        gofactor=gofactor-4;
        stopfactor=stopfactor+3;
    end
    if counter==10
        gofactor=gofactor-4;
        stopfactor=stopfactor+3;
    end
    
    if counter==14
        gofactor=gofactor-4;
        stopfactor=stopfactor+3;
    end
    
    timeint=gofactor*rand();

    tic
    while toc<=timeint
          pause(0.5)
          fprintf('GO')
    end
    fprintf('***STOP*** \n')


    timeint=stopfactor*rand();
    
    skiprnd=rand();
    if (1-skiprnd < skipfactor)
        %disp(skiprnd)
        timeint=0;
    end

    tic
    while toc<=timeint
          pause(0.5)
          fprintf('STOP')
    end
    fprintf('***GO*** \n')
    
    
end
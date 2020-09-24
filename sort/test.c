void test()
{
    
    Double_t a = -1.021390374;
    Double_t b = 5591.005348;
    Double_t c =  0;
    Double_t y = 0;
    Double_t counter = 0;
    Double_t x1, x2, y1, y2, m;
    x1=393;
    y1=1526;
    x2=468;
    y2=1255;
    m=(y1-y2)/(x1-x2);
    c=y1-m*x1;

    //Double_t beta = 0.019045;
    //Double_t gamma_energy = 4445;
    //Double_t angle = 180.0*3.142/180.0;
    //y = gamma_energy*(1.0 - beta*cos(angle))/(sqrt(1-beta*beta));
    for (Int_t i = 428; i<466; i++) {
        //y=a*i*i+b*i+c;
        y=m*i+c;
        counter=counter+y;
        
    }
    cout<<counter<<endl;
}




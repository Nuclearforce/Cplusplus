#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <cstdio>
#include "TSystem.h"

//start 'global'
TRandom3 *r = new TRandom3();
TFile f_hist("../histogram/14C/01_17_2018_7.root", "recreate");
//TFile f_hist("../histogram/13C_uncalibrated_before_experiment.root", "recreate");

TChain DataChain("data");
TTreeReader myReader(&DataChain);
TTreeReaderArray<unsigned int> TRenergy(myReader, "energy");
TTreeReaderArray<unsigned int> TRtime(myReader, "time");
TTreeReaderArray<unsigned int> TRmodchan(myReader, "modchan");

TChain Background("data");
TTreeReader myReader_background(&Background);
TTreeReaderArray<unsigned int> background_energy(myReader_background, "energy");
TTreeReaderArray<unsigned int> background_time(myReader_background, "time");
TTreeReaderArray<unsigned int> background_modchan(myReader_background, "modchan");

//extra
Bool_t for_testing = false;

//gate
Bool_t gate_on_excitation_spectra;
Bool_t particle_gamma;
Bool_t energy_sharing_ring;
Bool_t energy_sharing_sector;
Bool_t ring_sector_time_gate;
Double_t ring_time;
Double_t sector_time;
Double_t gamma_time;
Bool_t bad_gamma;


//lower and upper limit for excitation spectrum gate
Double_t lower_lim = -1;
Double_t upper_lim = 99999999;

//histograms
TH1F* hPixieADC[4][16];
//TH1F* hPixieADC_eff[4][16];
//TH1F* time_background_eff[4][16];
TH1F* hPixieADC_ab[4][16];
TH1F* hPixieADC_subtracted[4][16];
TH1F* hPixieADC_pg[4][16];
TH1F* backgroundADC[4][16];
TH1F* time_background_not_dopler[4][16];
TH1F* time_Ge;
TH1F* time_Cl;
TH1F* time_La;
TH1F* time_La_army;
TH1F* pg_time;
TH1F* ring_sector_time;
TH1F* time_background[4][16];
TH1F* time_Si;
TH1F* Clover1;
TH1F* Clover3;
TH1F* Clover1ab;
TH1F* Clover3ab;

//calibration
Double_t calib_a[4][16];
Double_t calib_b[4][16];
Double_t calib_c[4][16];
Double_t calib_d[4][16];
Double_t calib_e[4][16];

//Doppler correction
Double_t beta = 0.0;
Double_t theta[4][16];
Double_t phi[4][16];
Int_t ring_mod;
Int_t ring_chn;
Int_t sector_mod;
Int_t sector_chn;

//beta calculation
Double_t beta_ring_en=0.0;

//energy_sharing
Double_t energy_sharing_ring_en=0;
Double_t energy_sharing_ring_chn_use=0;
Double_t energy_sharing_ring_mod_use=0;
Double_t energy_sharing_ring_chn_ignore=0;
Double_t energy_sharing_ring_mod_ignore=0;
Double_t energy_sharing_sector_en=0;
Double_t energy_sharing_sector_chn_use=0;
Double_t energy_sharing_sector_chn_ignore=0;
Bool_t bad_particle;

//Add-back
Double_t ab_use_energy_cl1;
Double_t ab_use_crystal_cl1;
Double_t ab_ignore_crystal_cl1;
Bool_t ab_event_cl1=false;
Double_t ab_use_energy_cl3;
Double_t ab_use_crystal_cl3;
Double_t ab_ignore_crystal_cl3;
Bool_t ab_event_cl3=false;
Double_t test_counter=0;

//function declaration
void init_histos();
void read_calibrate_file();
void make_chain();
void first_loop();
void second_loop();
Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel);
Double_t do_doppler(Double_t en,Int_t gamma_mod, Int_t gamma_chn);
Double_t do_beta(Double_t ring_en);

void three_leaves_14()
{
    
    make_chain();//add trees to chain
    init_histos(); //initialise histograms
    read_calibrate_file(); //read calibration file
    
    // start timer for progress bar
    std::clock_t start;
    Double_t duration;
    Long64_t timer = 0;
    Long64_t max_event = DataChain.GetEntries();
    
    Bool_t stop_now = false;
    Long64_t max_event_testing = 2e6;
    start = std::clock();
    
    stop_now = false;
    
    ////////////////////// Loop over all entries of the TChain.////////////////////////////
    while (myReader.Next()) {
        first_loop();
        //cout<<"one"<<endl;
        //gate_on_excitation_spectra=true; //set to true if no excitation gate to be used
        //bad_particle=false;//set to false if no Si data
        //bad_gamma=false;//set to false if no Si data
        if(bad_particle==true||bad_gamma==true){
            //skip event
        }
        else if (gate_on_excitation_spectra==true){
            //cout<<"good event" <<endl;
            //cout<<"spam"<<endl;
            second_loop();
            //cout<<"two"<<endl;
        }
        
        
        if (for_testing==true){
            if(timer%100000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event_testing << " [%] \r";
                cout.flush();
            }
            timer++;
            if(timer==max_event_testing){
                stop_now=true;
            }
        }
        else{
            if(timer%50000000==0){
                cout << "Event #: " << timer  << " Progress: " << 100.0*timer/max_event << " [%] \r";
                cout.flush();
            }
            timer++;
        }
        if (stop_now==true){
            break;
        }
    }
    
    ////////////////write and close histogram file////////////////////
    Int_t j=0;
    Int_t i=0;
    for(j=0;j<4;j++){
        for(i=0;i<16;i++){
            if (j==2&&(i==0||i==1||i==2||i==3||i==4||i==5)){
                hPixieADC_pg[j][i]->Add(time_background[j][i],-(12.0/600.0));
                hPixieADC_ab[j][i]->Add(time_background_not_dopler[j][i],-(12.0/600.0));
                //hPixieADC_eff[j][i]->Add(time_background_eff[j][i],-(16.0/600.0));
            }
            else{
                hPixieADC_pg[j][i]->Add(time_background[j][i],-(12.0/600.0));
                hPixieADC_ab[j][i]->Add(time_background_not_dopler[j][i],-(12.0/600.0));
                //hPixieADC_eff[j][i]->Add(time_background_eff[j][i],-(16.0/600.0));
            }
        }
    }
    Clover1->Add(hPixieADC_pg[1][0],1);
    Clover1->Add(hPixieADC_pg[1][1],1);
    Clover1->Add(hPixieADC_pg[1][2],1);
    Clover1->Add(hPixieADC_pg[1][3],1);
    Clover1->Add(time_background[1][0],-(11.0/600.0));
    Clover1->Add(time_background[1][1],-(11.0/600.0));
    Clover1->Add(time_background[1][2],-(11.0/600.0));
    Clover1->Add(time_background[1][3],-(11.0/600.0));
    Clover3->Add(hPixieADC_pg[3][0],1);
    Clover3->Add(hPixieADC_pg[3][1],1);
    Clover3->Add(hPixieADC_pg[3][2],1);
    Clover3->Add(hPixieADC_pg[3][3],1);
    Clover3->Add(time_background[3][0],-(11.0/600.0));
    Clover3->Add(time_background[3][1],-(11.0/600.0));
    Clover3->Add(time_background[3][2],-(11.0/600.0));
    Clover3->Add(time_background[3][3],-(11.0/600.0));
    
    Clover1ab->Add(hPixieADC_ab[1][0],1);
    Clover1ab->Add(hPixieADC_ab[1][1],1);
    Clover1ab->Add(hPixieADC_ab[1][2],1);
    Clover1ab->Add(hPixieADC_ab[1][3],1);
    Clover1ab->Add(time_background_not_dopler[1][0],-(11.0/600.0));
    Clover1ab->Add(time_background_not_dopler[1][1],-(11.0/600.0));
    Clover1ab->Add(time_background_not_dopler[1][2],-(11.0/600.0));
    Clover1ab->Add(time_background_not_dopler[1][3],-(11.0/600.0));
    Clover3ab->Add(hPixieADC_ab[3][0],1);
    Clover3ab->Add(hPixieADC_ab[3][1],1);
    Clover3ab->Add(hPixieADC_ab[3][2],1);
    Clover3ab->Add(hPixieADC_ab[3][3],1);
    Clover3ab->Add(time_background_not_dopler[3][0],-(11.0/600.0));
    Clover3ab->Add(time_background_not_dopler[3][1],-(11.0/600.0));
    Clover3ab->Add(time_background_not_dopler[3][2],-(11.0/600.0));
    Clover3ab->Add(time_background_not_dopler[3][3],-(11.0/600.0));
    
    //cout<<gamma_amount_of_detectors_count<<endl;
    f_hist.Write();
    duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "Timer: " << duration << " [s] to analyze " << max_event << " events." << endl;
    f_hist.Close();
    
}

////////////////////////first loop, checks for add-back and gates//////////////////
void first_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
    Int_t sector_counter = 0;
    Int_t gamma_counter=0;
    Int_t ring_counter=0;
    Double_t time;
    Double_t en;
    Double_t calib_en;
    
    //gates
    gate_on_excitation_spectra = false;
    Double_t gamma_amount_of_detectors_check =0;
    ring_time=999999;
    sector_time=99999;
    gamma_time=99999;
    bad_gamma=false;
    
    //energy_sharing
    Double_t old_sharing_ring_en=0;
    Double_t old_sharing_ring_chn=0;
    Double_t old_sharing_ring_mod=0;
    Double_t old_sharing_sector_en=0;
    Double_t old_sharing_sector_chn=0;
    bad_particle=false;
    energy_sharing_ring=false;
    energy_sharing_sector=false;
    ring_sector_time_gate=false;
    
    //Add-back for clover 1
    Double_t ab_oldtime_cl1=0;
    Double_t ab_olden_cl1=0;
    Double_t ab_chn_cl1=0;
    Double_t count_cl1=0;
    Double_t ab_oldlocation_cl1=0;
    ab_event_cl1=false;
    
    //Add-back for clover 3
    Double_t ab_oldtime_cl3=0;
    Double_t ab_olden_cl3=0;
    Double_t ab_chn_cl3=0;
    Double_t count_cl3=0;
    Double_t ab_oldlocation_cl3=0;
    ab_event_cl3=false;
    
    for (entry_id = 0, entry_size = TRenergy.GetSize(); entry_id < entry_size; ++entry_id) {
        
        Int_t modchan=TRmodchan[entry_id];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time=TRtime[entry_id];
        en=TRenergy[entry_id];
        
        if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)){
            //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
            calib_en=en;
        }
        else{
            calib_en=do_calibrate(en,mod,chn);
        }
        //////////////////////////////////////////gates///////////////////////////////////////
        
        ///////////////energy sharing////////////////////////
        
        if( (mod==1 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==3 && (chn==7||chn==8||chn==9||chn==10||chn==12||chn==13||chn==14||chn==15)) ) {
            
            beta_ring_en=calib_en;
            
            
            if (ring_counter==0){
                ring_time=time;
                ring_mod = mod;
                ring_chn = chn;
                old_sharing_ring_en=calib_en;
                old_sharing_ring_chn=chn;
                old_sharing_ring_mod=mod;
            }
            
            if(ring_counter==1){
                //cout <<old_sharing_ring_chn-chn<<" " << old_sharing_ring_chn<< " " << chn << endl;
                if(abs(old_sharing_ring_chn-chn)==1){
                    energy_sharing_ring_en=old_sharing_sector_en+en;
                    energy_sharing_ring_chn_use=chn;
                    energy_sharing_ring_chn_ignore=old_sharing_ring_chn;
                    energy_sharing_ring_mod_use=mod;
                    energy_sharing_ring_mod_ignore=old_sharing_ring_mod;
                    energy_sharing_ring=true;
                    ring_mod = mod;
                    ring_chn = chn;
                    //cout <<old_sharing_ring_chn-chn<<" " << old_sharing_ring_chn<< " " << chn << endl;
                }
                else{
                    bad_particle=true;
                    //cout <<"bad_particle"<<endl;
                }
            }
            if(ring_counter>1||calib_en<10){
                bad_particle=true;
                //cout <<"bad_particle"<<endl;
            }
            
            /*if((mod==1&&(chn==8||chn==9))||(mod==2&&chn==8)){ //added to exclude rings during efficiency calculations
             bad_particle=true;
             }*/
            
            
            //////////excitation spectra gate (particle spectra gates)//////////////////
            if( (calib_en>lower_lim && calib_en<upper_lim)||(energy_sharing_ring_en>lower_lim && energy_sharing_ring_en<upper_lim &&energy_sharing_ring==true)){
                gate_on_excitation_spectra=true;
            }
            if((energy_sharing_ring_en<lower_lim||energy_sharing_ring_en>upper_lim)&&(energy_sharing_ring==true)){
                bad_particle=true;
            }
            
            
            ring_counter++;
            
        }
        
        if(mod==0){
            if (sector_counter==0){
                sector_time=time;
                old_sharing_sector_en=calib_en;
                old_sharing_sector_chn=chn;
                sector_mod = mod;
                sector_chn = chn;
            }
            if (sector_counter==1){
                if(abs(old_sharing_sector_chn-chn)==1){
                    energy_sharing_sector_en=old_sharing_sector_en+en;
                    energy_sharing_sector_chn_use=chn;
                    energy_sharing_sector_chn_ignore=old_sharing_sector_chn;
                    energy_sharing_sector=true;
                    sector_mod = mod;
                    sector_chn = chn;
                }
                else{
                    bad_particle=true;
                }
            }
            if(sector_counter>1 || calib_en < 10){
                bad_particle=true;
                
            }
            sector_counter++;
            
        }
        
        
        //////////////////////add-back////////////////////
        if (mod==1 && (chn==0 || chn==1 || chn==2 || chn==3) && calib_en < 10000 && calib_en>20){
            
            if ((time-ab_oldtime_cl1)<101.0 && count_cl1==1){
                ab_event_cl1=true;
                ab_use_energy_cl1=ab_olden_cl1+calib_en;
                ab_use_crystal_cl1=ab_oldlocation_cl1;
                ab_ignore_crystal_cl1=chn;
                
            }
            if(count_cl1>2||ab_olden_cl1+calib_en>10000){
                ab_event_cl1=false;
            }
            ab_oldtime_cl1=time;
            ab_olden_cl1=calib_en;
            ab_oldlocation_cl1=chn;
            count_cl1++;
        }
        
        if (mod==3 && (chn==0 || chn==1 || chn==2 || chn==3) && calib_en < 10000 && calib_en>20){
            
            if ((time-ab_oldtime_cl3)<101.0 && count_cl3==1){
                ab_event_cl3=true;
                ab_use_energy_cl3=ab_olden_cl3+calib_en;
                ab_use_crystal_cl3=ab_oldlocation_cl3;
                ab_ignore_crystal_cl3=chn;
            }
            if(count_cl3>2||ab_olden_cl3+calib_en>10000){
                ab_event_cl3=false;
            }
            ab_oldtime_cl3=time;
            ab_olden_cl3=calib_en;
            ab_oldlocation_cl3=chn;
            count_cl3++;
        }
        
        if( (mod==1 && (chn==0 || chn==1 || chn==2 || chn==3)) || (mod==3 && (chn==0 || chn==1 || chn==2 || chn==3||chn==5)) || (mod==2 && (chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)) ){
            if (gamma_counter==0){
                gamma_time=time;
            }
            if(ab_event_cl3==true||ab_event_cl1==true){
                bad_gamma=false;
            }
            else if(abs(gamma_amount_of_detectors_check-mod)>0 && gamma_counter>0){
                bad_gamma=true;
            }
            gamma_counter++;
            gamma_amount_of_detectors_check=mod;
        }
    }
    if( (ring_time-sector_time+200)>180 && (ring_time-sector_time+200)<220){
        ring_sector_time_gate=true;
    }
    
}


////////////////////////second loop fills histogram//////////////////
void second_loop(){
    
    Int_t entry_id;
    Int_t entry_size;
    Double_t time;
    Double_t en;
    Double_t calib_en;
    Double_t time_gate_lower=198;
    Double_t time_gate_upper=211;
    Double_t time_gate_lower_La=198;
    Double_t time_gate_upper_La=211;
    
    Double_t time_gate_si_lower=198;   //used for efficiency
    Double_t time_gate_si_upper=215;   //used for efficiency
    
    //Double_t time_gate_lower=0; //use when you don't want a time gate
    //Double_t time_gate_upper=1000;  //use when you don't want a time gate
    //ring_sector_time_gate=true; //set to true when you don't want a ring sector gate
    //ring_time=gamma_time;
    
    Double_t background_gate_lower=399;
    Double_t background_gate_upper=1001;
    Double_t fill_energy = 0.0;
    Double_t LaBr_uncalib_ds;
    Bool_t allready_filled;
    /*
     Double_t ab_use_energy_cl1;
     Double_t ab_use_crystal_cl1;
     Double_t ab_ignore_crystal_cl1;
     Double_t ab_use_energy_cl3;
     Double_t ab_use_crystal_cl3;
     Double_t ab_ignore_crystal_cl3;
     */
    
    
    if (ring_sector_time_gate==true){
        ring_sector_time->Fill(sector_time-ring_time+200);
    }
    for (entry_id = 0, entry_size = TRenergy.GetSize(); entry_id < entry_size; ++entry_id) {
        allready_filled = false;
        Int_t modchan=TRmodchan[entry_id];
        Int_t mod=modchan/16;
        Int_t chn=modchan%16;
        time=TRtime[entry_id];
        en=TRenergy[entry_id];
        if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)){
            //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
            //cout<< entry_id <<endl;
            calib_en=en;
        }
        else{
            calib_en=do_calibrate(en,mod,chn);
        }
        
        
        if( (mod==1 && (chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==2 && (chn==8||chn==9||chn==10||chn==11||chn==12||chn==13||chn==14||chn==15)) || (mod==3 && (chn==7||chn==8||chn==9||chn==10||chn==12||chn==13||chn==14||chn==15)) || mod==0 ){
            time_Si->Fill(time+200);
        }
        
        if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3||chn==5)) ){
            pg_time->Fill(time-ring_time+200);
        }
        
        if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3)) ){
            // if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_Cl->Fill(time+200);
            //}
        }
        
        if(mod==3&&chn==5){
            time_Ge->Fill(time+200);
        }
        
        if( (mod==1&&(chn==0||chn==1||chn==2||chn==3))||(mod==3&&(chn==0||chn==1||chn==2||chn==3||chn==5)) ){
            //if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_Ge->Fill(time+200);
            //}
        }
        if( mod==2&&(chn==0||chn==1||chn==2||chn==3) ){
            // if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_La->Fill(time+200);
            // }
        }
        
        if( mod==2&&(chn==4||chn==5) ){
            //  if ((time+200)>195 && (time+200)<210 && ring_sector_time_gate==true){
            time_La_army->Fill(time+200);
            // }
        }
        
        hPixieADC[mod][chn]->Fill(calib_en); //This is not callibrated LaBr due to doppler being done before calibration
        
        if (energy_sharing_sector==true && chn==energy_sharing_sector_chn_use && mod==0){
            
            /*if ((time+200)>time_gate_si_lower && (time+200)<time_gate_si_upper && ring_sector_time_gate==true){
             hPixieADC_eff[mod][chn]->Fill(energy_sharing_sector_en);
             }
             if ((time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true){
             time_background_eff[mod][chn]->Fill(energy_sharing_sector_en);
             }*/
            
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper && ring_sector_time_gate==true&&calib_en>6000){
                hPixieADC_ab[mod][chn]->Fill(energy_sharing_sector_en);
                hPixieADC_pg[mod][chn]->Fill(energy_sharing_sector_en);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true&&calib_en>6000){
                time_background[mod][chn]->Fill(energy_sharing_sector_en);
                time_background_not_dopler[mod][chn]->Fill(energy_sharing_sector_en);
            }
            allready_filled = true;
        }
        else if(energy_sharing_sector==true && chn==energy_sharing_sector_chn_ignore && mod==0){
            //do nothing
            allready_filled = true;
        }
        
        if (energy_sharing_ring==true && chn==energy_sharing_ring_chn_use && mod==energy_sharing_ring_mod_use){
            /*if ((time+200)>time_gate_si_lower && (time+200)<time_gate_si_upper && ring_sector_time_gate==true){
             hPixieADC_eff[mod][chn]->Fill(energy_sharing_ring_en);
             }
             if ((time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true){
             time_background_eff[mod][chn]->Fill(energy_sharing_ring_en);
             }*/
            
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper && ring_sector_time_gate==true&&calib_en>6000){
                hPixieADC_pg[mod][chn]->Fill(energy_sharing_ring_en);
                hPixieADC_ab[mod][chn]->Fill(energy_sharing_ring_en);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper && ring_sector_time_gate==true&&calib_en>6000){
                time_background[mod][chn]->Fill(energy_sharing_ring_en);
                time_background_not_dopler[mod][chn]->Fill(energy_sharing_ring_en);
            }
            allready_filled = true;
        }
        else if(energy_sharing_ring==true && chn==energy_sharing_ring_chn_ignore && mod==energy_sharing_ring_mod_ignore){
            //do nothing
            //cout<<"spam2"<<endl;
            allready_filled = true;
        }
        if (ab_event_cl1==true && chn==ab_use_crystal_cl1 && mod==1 && ring_sector_time_gate==true){
            beta=do_beta(beta_ring_en);
            fill_energy=do_doppler(ab_use_energy_cl1, mod, chn);
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper&&calib_en<10000){
                hPixieADC_pg[mod][chn]->Fill(fill_energy);
                hPixieADC_ab[mod][chn]->Fill(ab_use_energy_cl1);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper&&calib_en<10000){
                time_background[mod][chn]->Fill(fill_energy);
                time_background_not_dopler[mod][chn]->Fill(ab_use_energy_cl1);
            }
            
            allready_filled = true;
        }
        else if(ab_event_cl1==true && chn==ab_ignore_crystal_cl1 && mod==1){
            //do nothing
            allready_filled = true;
        }
        if (ab_event_cl3==true && chn==ab_use_crystal_cl3 && mod==3 && ring_sector_time_gate==true){
            beta=do_beta(beta_ring_en);
            fill_energy=do_doppler(ab_use_energy_cl3, mod, chn);
            
            if (abs(gamma_time-ring_time)<15 && (time+200)>time_gate_lower && (time+200)<time_gate_upper&&calib_en<10000 ){
                hPixieADC_pg[mod][chn]->Fill(fill_energy);
                hPixieADC_ab[mod][chn]->Fill(ab_use_energy_cl3);
            }
            if (abs(gamma_time-ring_time)<15 && (time+200)>background_gate_lower && (time+200)<background_gate_upper&&calib_en<10000){
                time_background[mod][chn]->Fill(fill_energy);
                time_background_not_dopler[mod][chn]->Fill(ab_use_energy_cl3);
            }
            allready_filled = true;
        }
        else if(ab_event_cl3==true && chn==ab_ignore_crystal_cl3 && mod==3){
            //do nothing
            allready_filled = true;
        }
        
        if(allready_filled==false && calib_en>20 && ring_sector_time_gate==true){
            if (energy_sharing_ring==true){
                ring_mod = energy_sharing_ring_mod_use;
                ring_chn = energy_sharing_ring_chn_use;
            }
            if (energy_sharing_sector==true){
                sector_chn=energy_sharing_sector_chn_use;
                sector_mod=0;
            }
            
            /* if ((time+200)>time_gate_si_lower && (time+200)<time_gate_si_upper){
             
             hPixieADC_eff[mod][chn]->Fill(calib_en);
             
             }
             if ((time+200)>background_gate_lower && (time+200)<background_gate_upper){
             time_background_eff[mod][chn]->Fill(calib_en);
             }*/
            
            if(mod == energy_sharing_ring_mod_ignore && chn == energy_sharing_ring_chn_ignore){
                //ignore
            }
            else if (abs(gamma_time-ring_time)<15){
                
                if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)&&calib_en<10000 && (time+200)> time_gate_lower_La && (time+200)<time_gate_upper_La){
                    //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
                    beta=do_beta(beta_ring_en);
                    LaBr_uncalib_ds=do_doppler(calib_en, mod, chn);
                    fill_energy=do_calibrate(LaBr_uncalib_ds,mod,chn);
                    calib_en=do_calibrate(calib_en,mod,chn);
                    hPixieADC_pg[mod][chn]->Fill(fill_energy);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                }
                else if(((mod==1&&(chn==0 || chn==1 || chn==2 || chn==3))||(mod==3&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==5)))&& calib_en<10000   && (time+200)>time_gate_lower && (time+200)<time_gate_upper){
                    beta=do_beta(beta_ring_en);
                    fill_energy=do_doppler(calib_en, mod, chn);
                    hPixieADC_pg[mod][chn]->Fill(fill_energy);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                }
                else if((time+200)>time_gate_lower && (time+200)<time_gate_upper&&calib_en>6000){ //i.e. if it is the silicon
                    //fill_energy=calib_en;
                    hPixieADC_pg[mod][chn]->Fill(calib_en);
                    hPixieADC_ab[mod][chn]->Fill(calib_en);
                }
                if(mod==2&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==4 || chn==5)&&calib_en<10000 && (time+200)>background_gate_lower && (time+200)<background_gate_upper){
                    //for the LaBr we want to doppler shift first and then calibrate. For the Ge we want to calibrate, Add-back then doppler shift.
                    beta=do_beta(beta_ring_en);
                    LaBr_uncalib_ds=do_doppler(calib_en, mod, chn);
                    fill_energy=do_calibrate(LaBr_uncalib_ds,mod,chn);
                    calib_en=do_calibrate(calib_en,mod,chn);
                    time_background[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler[mod][chn]->Fill(calib_en);
                }
                else if(((mod==1&&(chn==0 || chn==1 || chn==2 || chn==3))||(mod==3&&(chn==0 || chn==1 || chn==2 || chn==3 || chn==5)))&&calib_en<10000 && (time+200)>background_gate_lower && (time+200)<background_gate_upper){
                    beta=do_beta(beta_ring_en);
                    fill_energy=do_doppler(calib_en, mod, chn);
                    time_background[mod][chn]->Fill(fill_energy);
                    time_background_not_dopler[mod][chn]->Fill(calib_en);
                    
                }
                else if((time+200)>background_gate_lower && (time+200)<background_gate_upper&&calib_en>6000){ //i.e. if it is the silicon
                    //fill_energy=calib_en;
                    time_background[mod][chn]->Fill(calib_en);
                    time_background_not_dopler[mod][chn]->Fill(calib_en);
                }
                
            }
            
        }
        
    }
}

void read_calibrate_file() {
    Double_t a = 0.0;
    Double_t b = 1.0;
    Double_t c = 1.0;
    Double_t d = 0.0;
    Double_t e = 0.0;
    Double_t f = 0.0;
    Double_t g = 0.0;
    FILE *calib = fopen("Calibration_before_exp_14C.txt","r");
    //FILE *calib = fopen("Calibration_empty.txt","r");
    for(int m=0;m<4;m++){
        for(int n=0;n<16;n++){ //The format is: offset(a) slope^1(b) slope(c)^2 phi(d) theta(e)
            fscanf(calib, "%lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &f, &g, &d, &e);
            calib_a[m][n]=a;
            calib_b[m][n]=b;
            calib_c[m][n]=c;
            calib_d[m][n]=f;
            calib_e[m][n]=g;
            phi[m][n]=d;
            theta[m][n]=e;
        }
    }
    
}

Double_t do_calibrate(Double_t raw_energy, Int_t module, Int_t channel){
    Double_t random_num;
    Double_t cal_energy;
    Double_t offset = calib_a[module][channel];
    Double_t slope1 = calib_b[module][channel];
    Double_t slope2 = calib_c[module][channel];
    Double_t slope3 = calib_d[module][channel];
    Double_t slope4 = calib_e[module][channel];
    random_num=r->Uniform(-0.5,0.5);
    cal_energy=(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*slope4+(raw_energy+random_num)*(raw_energy+random_num)*(raw_energy+random_num)*slope3+(raw_energy+random_num)*(raw_energy+random_num)*slope2 + (raw_energy+random_num)*slope1 + offset;
    //cal_energy=raw_energy+random_num;
    //if (module==2&&channel==14){
    //    cout << raw_energy << " " << cal_energy << endl;
    //}
    return cal_energy;
}


Double_t do_beta(Double_t ring_en){
    Double_t beta;
    beta=sqrt(((ring_en/1000.0)*2.0)/(931.5*14.0));
    /*if (beta==0){
     cout<<ring_en<<endl;
     }*/
    //cout<<beta<<endl;
    return beta;
}


Double_t do_doppler(Double_t en, Int_t gamma_mod, Int_t gamma_chn){
    Double_t dop_cor_en;
    Double_t cos_alpha;
    //cos_alpha=sin(theta[ring_mod][ring_chn])*sin(theta[gamma_mod][gamma_chn])*cos(phi[sector_mod][sector_chn]-phi[gamma_mod][gamma_chn]) + cos(theta[ring_mod][ring_chn])*cos(theta[gamma_mod][gamma_chn]);
    cos_alpha=cos(theta[gamma_mod][gamma_chn]-theta[ring_mod][ring_chn]);
    dop_cor_en=(en)*(1.0+beta*cos_alpha)/(sqrt(1-beta*beta));
    if (gamma_mod==3&&gamma_chn==5){
        //cout << "cos alpha " << cos_alpha << " energy before " << en << " doppler shifted energy " << dop_cor_en << " theta_p " << theta[ring_mod][ring_chn]  << " phi_p " << phi[sector_mod][sector_chn] << " phi_gamma " << phi[gamma_mod][gamma_chn]<< " theta_gamma " << theta[gamma_mod][gamma_chn] << " beta " << beta << endl;
        
    }
    return dop_cor_en;
}


void init_histos() {
    //calibrated
    /*Int_t gamma_bins=10000;
     Int_t gamma_min=0;
     Int_t gamma_max=9999;
     Int_t particle_bins=36000;
     Int_t particle_min=0;
     Int_t particle_max=35999;*/
    
    //uncalibrated
    Int_t gamma_bins=32768;
    Int_t gamma_min=0;
    Int_t gamma_max=32768;
    Int_t particle_bins=32768;
    Int_t particle_min=0;
    Int_t particle_max=32768;
    
    
    time_Ge = new TH1F("time_Ge_single","",1500,0,1500);
    time_Cl = new TH1F("time_Cl","",1500,0,1500);
    time_La = new TH1F("time_La","",1500,0,1500);
    time_La_army = new TH1F("time_La_army","",1500,0,1500);
    ring_sector_time = new TH1F("time_ring_sector","",1500,0,1500);
    pg_time = new TH1F("time_particle_gamma","",1500,0,1500);
    time_Si = new TH1F("time_Si","",1500,0,1500);
    Clover1 = new TH1F("Clover_1","",particle_bins,particle_min,particle_max);
    Clover3 = new TH1F("Clover_3","",particle_bins,particle_min,particle_max);
    Clover1ab = new TH1F("Clover_1ab","",particle_bins,particle_min,particle_max);
    Clover3ab = new TH1F("Clover_3ab","",particle_bins,particle_min,particle_max);
    
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC_pg[j][i] = new TH1F(Form("ADC_complete_%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC_ab[j][i] = new TH1F(Form("ADC_ab%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            hPixieADC[j][i] = new TH1F(Form("ADC%d_ch%02d",j,i),"",50000,0,49999);
        }
    }
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            time_background[j][i] = new TH1F(Form("time_background%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    
    for(Int_t j=0;j<4;j++){
        
        for(Int_t i=0;i<16;i++){
            time_background_not_dopler[j][i] = new TH1F(Form("time_background_not_dopler%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
        }
    }
    
    /*for(Int_t j=0;j<4;j++){
     
     for(Int_t i=0;i<16;i++){
     hPixieADC_eff[j][i] = new TH1F(Form("ADC_eff_%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
     }
     }
     
     for(Int_t j=0;j<4;j++){
     
     for(Int_t i=0;i<16;i++){
     time_background_eff[j][i] = new TH1F(Form("time_background_eff%d_ch%02d",j,i),"",particle_bins,particle_min,particle_max);
     }
     }*/
    
    
}

void make_chain(){
    
    //////14C/////////////
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0675*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0676*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0677*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0678*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0679*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-068*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-069*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-070*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-071*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-072*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-073*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-074*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-075*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-076*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-077*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-078*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0790*.root");
    DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0791*.root");
    /* DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0809*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0810*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0811*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0812*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0813*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0814*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0815*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0816*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0817*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0818*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0822*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0823*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0824*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0825*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0826*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0827*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0828*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0829*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-083*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-084*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-085*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-086*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-087*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-088*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-089*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0901*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0902*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0903*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0904*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0905*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0906*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0907*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0908*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0910*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0911*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0912*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0913*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0914*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0915*.root");
     DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0916*.root");*/
    
    //DataChain.Add("/Volumes/CBRITS4TB/c0dea36f-31d2-44d7-8418-e1a6319e7468/14C/run-0689-02*.root");   //one for testing (900mb)
    
}


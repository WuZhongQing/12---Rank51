#include<iostream>
#include<vector>
#include<math.h>
#include<float.h>
#include<fstream>

using namespace std;

double v=5, H=10, Dinter=80, Dintra=90,d=125,D=70;
double tf=0.1;

vector<double> compute_ABC(vector<double> &start_site, vector<double> &end_site){
    /*compute the line AX+BY+C=0*/
   double A=1.0,B,C;
   B = -(start_site[1]-end_site[1])/(start_site[0]-end_site[0]);
   C = -start_site[1] - start_site[0]*B;
   return {A,B,C};
}

void show(vector<double> &a){
    for(auto x:a){
        cout<<x<<" ";
    }
    cout<<endl;
}
double compute_distance_site(vector<double> &site_A, vector<double> &site_B){
    /*compute the distance between each site*/
    
    double d = pow(site_A[0]-site_B[0],2) + pow(site_A[1]-site_B[1],2) + pow(site_A[2]-site_B[2],2);
    return pow(d, 1.0/2);
}



vector<vector<int>> find_can_send(vector<double> &center_site, double distance, double t, double Dh, vector<int> &m_n){
    vector<vector<int>> res;
    double time = tf + distance/10000.0 + t;

    double d_cir = sqrt(distance*distance - Dh*Dh); //compute the radius of the cirle
    double xmax = center_site[0] + d_cir, xmin = center_site[0] - d_cir;
    double ymax = center_site[1] + d_cir, ymin = center_site[1] - d_cir;
    int nmax = (int)(ymax/Dinter), nmin = (int)(ymin/Dintra);
    int mmax = (int)((xmax-v*time)/Dintra),mmin = (int)((xmin-v*time)/Dintra);

    if(nmin<0){
        nmin-=1;
    } 
    if(mmax-mmin==0 && nmax-mmin==0){
        cout<<"OMG"<<endl;
    }
    for(int m=mmin;m<=mmax;m++){
        for(int n=nmin;n<=nmax;n++){
            double x = v*time+m*Dintra, y = n*Dinter;
            vector<double> site_temp = {x,y,H};
            // show(center_site);
            // cout<<m_n[0]<<" "<<m_n[1]<<endl;
            if(compute_distance_site(site_temp, center_site) <= distance && (m!=m_n[0] || n!=m_n[1])){
                res.push_back({m,n});
            }
        }
    }

    return res;
}



double compute_dis_point_line(vector<double> &point, vector<double> &ABC){
    return abs(-(point[1]+ABC[2])/ABC[1] - point[0] );
    // return abs(ABC[0]*point[0]+ABC[1]*point[1]+ABC[2])/sqrt(pow(ABC[0],2)+pow(ABC[1],2));
}

vector<double> close_to_end(vector<double> &start_site, vector<vector<int>> &can_send_end, vector<double> &end_site, double start_time, double end_time, vector<int> &m_n){
    vector<double> res;
    int m_finally,n_finally;
    double time_temp = (end_time-start_time)/1000.0;
    double s=FLT_MAX;
    for(int i=1;i<=1000;i++){
        for(auto site:can_send_end){
            // cout<<"work"<<endl;
            int m=site[0],n=site[1];
            double x = v*(start_time + time_temp*i)+m*Dintra, y = n*Dinter, h=H;
            vector<double> temp = {x,y,h};
            double dis = compute_distance_site(temp, end_site);
            double dis_s_end =  compute_distance_site(start_site, temp);
            if(dis<=D){
                if(dis+ dis_s_end <s){
                    s = dis+ dis_s_end;
                    // res = site;
                    res = temp;
                    m_finally = m;
                    n_finally = n;
                }
            }
        }
    }
    m_n[0] = m_finally;
    m_n[1] = n_finally;
    return res;

}

vector<double> get_next_v1(vector<double> &ABC, vector<double> &start_site, vector<double> &end_site, vector<vector<int>> &can_send,double start_time, double end_time, vector<int> &m_n){
    vector<vector<int>> can_send_end;
    double min_d = FLT_MAX;
    vector<double> res;
    int m_end,n_end;

    double time_temp = (end_time-start_time)/100.0;
    for(int i=0;i<=100;i++){
        for(auto site:can_send){
            int m=site[0],n=site[1];
            double x = v*(start_time   + time_temp*i)+m*Dintra, y = n*Dinter,h=H;
            vector<double> temp = {x,y,h};
            double dis = compute_distance_site(temp, end_site);
            double dis2 = compute_distance_site(temp, start_site);
            if(dis<=D){
                can_send_end.push_back({m,n});
                // res = temp;
                // m_n[0] = m;
                // m_n[1] = n;
                // return res;
            }
            if(dis <= min_d && dis2<d){
                min_d = dis;
                res = temp;
                m_end = m;
                n_end = n;
            }

        }
    }
    m_n[0] = m_end;
    m_n[1] = n_end;
    if(can_send_end.size()!=0){
        vector<double> m_n_temp =  close_to_end(start_site, can_send_end, end_site, start_time, end_time, m_n);
        res = m_n_temp;


    }
    return res;
}

vector<double> get_next_v2(vector<double> &ABC, vector<double> &start_site,vector<double> &end_site, vector<vector<int>> &can_send,double start_time, double end_time, vector<int> &m_n){
    vector<vector<int>> can_send_end;
    double min_dis = FLT_MAX;
    vector<double> res;

    int m_finally,n_finally;
    int m_finally_close, n_finally_close;
    double d_start_end = compute_distance_site(start_site, end_site);
    double time_step= (end_time-start_time)/100.0;
    bool flag=true;
    if(can_send.size() <= 2){
        flag = false;
    }
    for(int i=1;i<=100;i++){
        for(auto site:can_send){
            int m=site[0],n=site[1];
            double x = v*(start_time + time_step*i)+m*Dintra, y = n*Dinter,h=H;
            vector<double> temp = {x,y,h};

            double dis = compute_distance_site(temp, end_site);
            double dis2 = compute_distance_site(temp, start_site);
            double dis_point = compute_dis_point_line(temp,ABC);
            if(dis<=D){
                can_send_end.push_back({m,n});
            //     res = temp;
            //     m_n[0] = m;
            //     m_n[1] = n;
            //     return res;
            }
            if(dis < d_start_end && dis_point < min_dis && flag && dis2<d){

                min_dis = dis_point;
                res = temp;
                m_finally = m;
                n_finally = n;

            }else if(dis < min_dis && !flag && dis2<d){
                min_dis = dis;
                res = temp;
                m_finally = m;
                n_finally = n;
            }

        }
    }
    // cout<<min_d<<endl;
    m_n[0] = m_finally;
    m_n[1] = n_finally;
    if(can_send_end.size()!=0){
        vector<double> m_n_temp =  close_to_end(start_site, can_send_end, end_site, start_time, end_time, m_n);
        res = m_n_temp;
    }
    return res;
}

float check_time_is_legal(float start_time, float end_time, vector<float> &start_site, vector<float> &end_site){
    if(end_time - start_time < tf + compute_distance_site(start_site, end_site)/10000.0){
        return start_time + tf + compute_distance_site(start_site, end_site)/10000.0 + 0.0001;
    }
    
    return end_time+0.0001;
}

int main(int argc, char const *argv[])
{

    vector<double> start_time = {0, 4.7, 16.4};
    vector<double> start_site = {45.73, 45.26, 0};
    vector<vector<double>> end_sites = {{1200, 700,0},{-940, 1100, 0}};
    ofstream ofs;
    ofs.open("./result.txt", ios::app);
    
    for(auto time:start_time){
        int index_end_1=1,index_start_1=0;
        int index_end_2 = 1, index_start_2=0;
        
        for(auto end_site:end_sites){
            vector<double> ABC = compute_ABC(start_site, end_site);// compute A B C


            string line1 = to_string(time)+","+to_string(index_start_1);
            string line2 = to_string(time)+","+to_string(index_start_2);;
            string line_each1, line_each2;
            double time_temp_1 = time, time_temp_2 = time;
            double time_1 = time, time_2=time;
            
            bool flag_1 = true;
            bool flag_2 =  true;
            vector<double> start_1 = start_site;
            vector<double> start_2 = start_site;
            double t1=0,t2=0;
            vector<int> m_n_1 = {INT8_MAX, INT8_MAX};
            vector<int> m_n_2 = {INT8_MAX, INT8_MAX};
            while(flag_1){
                double  Dh = start_1[2]==0 ? H:0;
                double dis = start_1[2]==0 ? D:d;
                
                vector<vector<int>> can_send = find_can_send(start_1, dis, time_temp_1, Dh, m_n_1);
                double start_time = time_temp_1;
                time_temp_1 += tf + dis/10000.0;
                vector<double> temp = start_1;
                start_1 = get_next_v2(ABC, start_1, end_site, can_send, start_time, time_temp_1, m_n_1);

                // cout<<start[0]<<" "<<start[1]<<" "<<start[2]<<endl;
                double S =  compute_distance_site(temp, start_1);
                time_1 += tf + S/10000.0;
                t1+=tf + S/10000.0;
                line_each1+="("+to_string(time_1)+","+to_string(m_n_1[0])+","+to_string(m_n_1[1])+"),";

                if(compute_distance_site(start_1, end_site)<=D){
                    flag_1=false;
                    t1+= tf + compute_distance_site(start_1, end_site)/10000.0;
                    line1+= ","+to_string(index_end_1++)+","+to_string(t1);
                    line_each1 = line_each1.substr(0, line_each1.size()-1);  
                }
            }
            //********************another **********************//
            while(flag_2){
                double  Dh = start_2[2]==0 ? H:0;
                double dis = start_2[2]==0 ? D:d;
                
                vector<vector<int>> can_send = find_can_send(start_2, dis, time_temp_2, Dh, m_n_2);
                double start_time = time_temp_2;
                time_temp_2 += tf + dis/10000.0;
                vector<double> temp = start_2;
                start_2 = get_next_v1(ABC, start_2, end_site, can_send, start_time, time_temp_2, m_n_2);

                // cout<<start[0]<<" "<<start[1]<<" "<<start[2]<<endl;
                double S =  compute_distance_site(temp, start_2);
                time_2 += tf + S/10000.0;
                t2+=tf + S/10000.0;
                line_each2+="("+to_string(time_2)+","+to_string(m_n_2[0])+","+to_string(m_n_2[1])+"),";

                if(compute_distance_site(start_2, end_site)<=D){
                    flag_2=false;
                    t2+= tf + compute_distance_site(start_2, end_site)/10000.0;
                    line2+= ","+to_string(index_end_2++)+","+to_string(t2);
                    line_each2 = line_each2.substr(0, line_each2.size()-1);  
                }
            }
            // cout<<"t1:"<<t1<<endl;
            // cout<<"t2:"<<t2<<endl;
            if(t1<t2){

                ofs<<line1<<endl;
                cout<<line1<<endl;
                ofs<<line_each1<<endl;
                cout<<line_each1<<endl;
            }else{

                ofs<<line2<<endl;
                cout<<line2<<endl;
                ofs<<line_each2<<endl;
                cout<<line_each2<<endl; 
            }

        }
    }
    ofs.close();
    return 0;
}

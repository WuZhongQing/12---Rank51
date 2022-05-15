#include<iostream>
#include<vector>
#include<math.h>
#include<float.h>
#include<fstream>
#include<algorithm>
#include<unordered_map>
#include<sstream>
#include<iomanip>

using namespace std;

float v=5, H=10, Dinter=80, Dintra=90, d=115, D=70;
float tf=0.1;
float sig_num = 10;


vector<vector<float>> sky_stations = {{-614,1059,24},{-934, 715,12},{1073,291,37},{715,129,35},{186,432,21},{-923,632,37},{833,187,24},{-63,363,11}};


int find_index(vector<vector<float>> &sky_stations, vector<float> &sky_station){
    return find(sky_stations.begin(), sky_stations.end(), sky_station) - sky_stations.begin();
}

vector<float> compute_ABC(vector<float> &start_site, vector<float> &end_site){
    /*compute the line AX+BY+C=0*/
   float A=1.0,B,C;
   B = -(start_site[1]-end_site[1])/(start_site[0]-end_site[0]);
   C = -start_site[1] - start_site[0]*B;
   return {A,B,C};
}

float compute_dis_point_line(vector<float> &point, vector<float> &ABC){
    return abs(-(point[1]+ABC[2])/ABC[1] - point[0]);
    // return abs(ABC[0]*point[0]+ABC[1]*point[1]+ABC[2])/sqrt(pow(ABC[0],2)+pow(ABC[1],2));
}

float compute_distance_site(vector<float> &site_A, vector<float> &site_B){
    /*compute the distance between each site*/
    
    float d = pow(site_A[0]-site_B[0],2) + pow(site_A[1]-site_B[1],2) + pow(site_A[2]-site_B[2],2);
    return pow(d, 1.0/2);
}

void show_site(vector<float> &start_site_){
    for(auto num:start_site_ ){
        cout<<num<<" ";
    }
    cout<<endl;
}


vector<vector<int>> find_can_send_v2(vector<float> &start_site,float time_, vector<int> &last_mn){
    vector<vector<int>> res;
    if(last_mn[0]<100 && last_mn[0]>-100){
        res.push_back({last_mn[0], last_mn[1]+1});
        res.push_back({last_mn[0], last_mn[1]-1});
        res.push_back({last_mn[0]+1, last_mn[1]});
        res.push_back({last_mn[0]-1, last_mn[1]});
    }else{
        float time = tf +d/10000.0+time_;

        float d_cir = sqrt(d*d - (start_site[2]-H)*(start_site[2]-H)); //compute the radius of the cirle
        float xmax = start_site[0] + d_cir, xmin = start_site[0] - d_cir;
        float ymax = start_site[1] + d_cir, ymin = start_site[1] - d_cir;
        int nmax = (int)(ymax/Dinter), nmin = (int)(ymin/Dintra);
        int mmax = (int)((xmax-v*time)/Dintra),mmin = (int)((xmin-v*time)/Dintra);

        if(nmin<0){
            nmin-=1;
        } 
        for(int m=mmin;m<=mmax;m++){
            for(int n=nmin;n<=nmax;n++){
                float x = v*time+m*Dintra, y = n*Dinter;
                vector<float> site_temp = {x,y,H};
                if(compute_distance_site(site_temp, start_site) <= d){
                    res.push_back({m,n});
                }
            }
        }
    }
    return res;
}


vector<vector<int>> find_first(float t, float distance, vector<float> &start_site,vector<float> &end_site){
    /*find who is the first step*/
    vector<vector<int>> res_temp;
    float time = tf + distance/10000.0 + t;
    float d_cir = sqrt(distance*distance - H*H); //compute the radius of the cirle
    float xmax = start_site[0] + d_cir, xmin = start_site[0] - d_cir;
    float ymax = start_site[1] + d_cir, ymin = start_site[1] - d_cir;
    int nmax = (int)(ymax/Dinter), nmin = (int)(ymin/Dintra);
    int mmax = (int)((xmax-v*time)/Dintra),mmin = (int)((xmin-v*time)/Dintra);
    for(int m=mmin;m<=mmax;m++){
        for(int n=nmin;n<=nmax;n++){
            float x = v*time+m*Dintra, y = n*Dinter;
            vector<float> site_temp = {x,y,H};
            if(compute_distance_site(site_temp, start_site) <= distance){
                res_temp.push_back({m,n});
            }
        }
    }
    vector<vector<int>> res;
    while(true){
        if(res_temp.size()==1){
            break;
        }
        float dis = FLT_MAX;
        vector<int> temp;
        for(auto mn:res_temp){
            vector<float> site_temp = {v*time + mn[0]*Dintra, mn[1]*Dinter, H};
            if(compute_distance_site(site_temp, end_site) <= dis){
                dis = compute_distance_site(site_temp, end_site);
                temp = mn;

            }

        }
        res_temp.erase(find(res_temp.begin(), res_temp.end(), temp));
        res.push_back(temp);
    }
    res.push_back(res_temp[0]);

    return res;
}

vector<float> close_to_end(vector<float> &start_site, vector<vector<int>> &can_send_end, vector<float> &end_site, float start_time, float end_time, vector<int> &m_n){
    vector<float> res;
    int m_finally,n_finally;
    float time_temp = (end_time-start_time)/1000.0;
    float s=FLT_MAX;
    for(int i=1;i<=1000;i++){
        for(auto site:can_send_end){
            // cout<<"work"<<endl;
            int m=site[0],n=site[1];
            float x = v*(start_time + time_temp*i)+m*Dintra, y = n*Dinter, h=H;
            vector<float> temp = {x,y,h};
            float dis = compute_distance_site(temp, end_site);
            float dis_s_end =  compute_distance_site(start_site, temp);
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



vector<vector<float>> find_sky_station_up(vector<float> &ABC, vector<float> &start_site,vector<float> &end_site){
    vector<vector<float>> res_temp;
    float max_x = max(start_site[0], end_site[0]);
    float min_x = min(start_site[0], end_site[0]);
    float max_y = max(start_site[1], end_site[1]);
    float min_y = min(start_site[1], end_site[1]);
    
    for(auto sky_station:sky_stations){
        if(sky_station[0]<max_x &&sky_station[0]>min_x && sky_station[1] < max_y && sky_station[1] > min_y){
            float y = -(ABC[1]*sky_station[0]+ABC[2]);
            if(y<sky_station[1]){
                res_temp.push_back(sky_station);
            }
        }
    }
    if(res_temp.size()<=1){
        return res_temp;
    }
    vector<vector<float>> res;
    while (res_temp.size()>1){
        
        float dis = FLT_MAX;
        vector<float> m;
        for(auto site:res_temp){
            if(compute_distance_site(site, end_site) < dis){
                dis = compute_distance_site(site, end_site);
                m = site;
            }
        }
        res.push_back(m);
        res_temp.erase(find(res_temp.begin(), res_temp.end(), m));
    }
    res.push_back(res_temp[0]);
    return res;
    
}
vector<vector<float>> find_sky_station_down(vector<float> &ABC, vector<float> &start_site,vector<float> &end_site){
    vector<vector<float>> res_temp;
    float max_x = max(start_site[0], end_site[0]);
    float min_x = min(start_site[0], end_site[0]);
    float max_y = max(start_site[1], end_site[1]);
    float min_y = min(start_site[1], end_site[1]);
    
    for(auto sky_station:sky_stations){
        if(sky_station[0]<max_x &&sky_station[0]>min_x && sky_station[1] < max_y && sky_station[1] > min_y){
            float y = -(ABC[1]*sky_station[0]+ABC[2]);
            if(y>sky_station[1]){
                res_temp.push_back(sky_station);
            }
        }
    }
    if(res_temp.size()<=1){
        return res_temp;
    }
    vector<vector<float>> res;
    while (res_temp.size()>1){
        float dis = FLT_MAX;
        vector<float> m;
        for(auto site:res_temp){
            if(compute_distance_site(site, end_site) < dis){
                dis = compute_distance_site(site, end_site);
                m = site;
            }
        } 
        res.push_back(m);
        res_temp.erase(find(res_temp.begin(), res_temp.end(), m));
    }
    res.push_back(res_temp[0]);
    // for(auto site:res){
    //     show_site(site);
    // }
    return res;
}


bool check_inside(vector<float> &start_site,vector<float> &sky_station, vector<float> &end_site){
    float x_max = max(start_site[0], end_site[0]);
    float x_min = min(start_site[0], end_site[0]);
    float y_max = max(start_site[1], end_site[1]);
    float y_min = min(start_site[1], end_site[1]);
    if(sky_station[0] < x_max-1 && sky_station[0] > x_min+1 && sky_station[1] < y_max-1 && sky_station[1] > y_min+1){
        return true;
        
    }
    return false;
}

vector<float> check_can_send_sky(vector<float> &start_site, vector<float> &sky_station, vector<int> &mn_){
    int index = find_index(sky_stations, sky_station);
    mn_[0] = 100+index;
    mn_[1] = 100+index;
    return sky_station;
}
vector<float> find_most_close(vector<float> &start_site, vector<vector<int>> &can_send,float time, vector<float> &end_site, vector<int> &mn_,vector<vector<int>> &can_send_end){

    float dis = FLT_MAX;
    vector<float> res;
    for(auto mn:can_send){
        vector<float> temp  = {v*time+mn[0]*Dintra, mn[1]*Dinter, H};
        float cost_time = tf + compute_distance_site(temp, start_site) / 10000.0;
        vector<float> res_temp = {v*(cost_time+time)+mn[0]*Dintra, mn[1]*Dinter, H};
        if(compute_distance_site(res_temp, end_site)<=D){
            can_send_end.push_back(mn);
        }

        if(compute_distance_site(res_temp, end_site) < dis){
            dis = compute_distance_site(res_temp, end_site);
            res = res_temp;
            mn_[0] = mn[0];
            mn_[1]  =mn[1];
        }
    }
    
    
    return res;
}
vector<float> find_most_close_v2(vector<float> &start_site, vector<vector<int>> &can_send,float time, vector<float> &end_site, vector<int> &mn_,vector<vector<int>> &can_send_end, unordered_map<string, vector<pair<float, float>>> &record_time){

    float dis = FLT_MAX;
    vector<float> res;
    for(auto mn:can_send){
        vector<float> temp  = {v*time+mn[0]*Dintra, mn[1]*Dinter, H};
        float cost_time = tf + compute_distance_site(temp, start_site) / 10000.0;
        vector<float> res_temp = {v*(cost_time+time)+mn[0]*Dintra, mn[1]*Dinter, H};

        string query = to_string((int)mn_[0])+to_string((int)mn_[1])+to_string((int)mn[0])+to_string((int)mn[1]);

        if(compute_distance_site(res_temp, end_site)<=D && record_time[query].size()==0){
            can_send_end.push_back(mn);
        }

        if(compute_distance_site(res_temp, end_site) < dis && record_time[query].size()==0){
            dis = compute_distance_site(res_temp, end_site);
            res = res_temp;
            mn_[0] = mn[0];
            mn_[1]  =mn[1];
        }
    }
    
    
    return res;
}

vector<float> get_site(vector<float> &start_site, vector<float> &end_site, string mode, vector<float> &ABC){
    //mode: up down
    //ABC y+Bx+C=0
    float x,y,h;
    if(mode=="up"){
        y = max(end_site[1], start_site[1]);
        if(-ABC[1]>0){
            x = min(end_site[0], start_site[0]);
        }else{
            x = max(end_site[0], start_site[0]);
        }
    }else{
        y = min(end_site[1], start_site[1]);
        if(-ABC[1]>0){
            x = max(end_site[0], start_site[0]);
        }else{
            x = min(end_site[0], start_site[0]);
        }
    }
    h = H;
    return {x,y,h};
}

bool is_good(vector<float> start_site, vector<float> sky, vector<float> &end_site, float start_time,vector<int> &m_n){
    vector<float> start_site_ = start_site;
    
    vector<int> mn = m_n;
    int index = find_index(sky_stations, sky);
    vector<int> sky_mn = {100 + index, 100+index};

    vector<vector<int>> can_send_end;
    vector<vector<int>> can_send = find_can_send_v2(sky, start_time, sky_mn);
    vector<float> most_close = find_most_close(sky, can_send, start_time, end_site, mn, can_send_end);

    float dis1 = compute_distance_site(start_site, sky) + compute_distance_site(most_close, sky) + 1000;
    float dis2 = 0;
    bool flag=true;
    while(flag){
        can_send = find_can_send_v2(start_site_, start_time, mn);
        vector<float> start_site_temp = start_site_;
        
        start_site_ = find_most_close(start_site_, can_send, start_time, most_close, mn, can_send_end);
        start_time += 0.1+85/10000.0;
        dis2 += compute_distance_site(start_site_, start_site_temp) + 1000;
        if(compute_distance_site(start_site_, most_close) <= d){
            flag = false;
            dis2+= compute_distance_site(start_site_, most_close) +1000;
        }
    }
    if(dis2>dis1){
        return true;
    }
    return false;
}

vector<float> get_next_v3(vector<vector<float>> &skys, vector<float> &start_site,vector<float> &end_site, vector<vector<int>> &can_send, float start_time, float end_time, vector<int> &m_n){
    vector<vector<int>> can_send_end;
    vector<float> res;
    for(int i=skys.size()-1;i>=0;i--){
        vector<float> sky_station = skys[i];
        if(check_inside(start_site, sky_station, end_site)){
            
            if(compute_distance_site(start_site, sky_station) <= d){
                bool flag = is_good(start_site,sky_station,end_site,start_time,m_n);
                vector<float> most_close;
                if(flag){
                    most_close = check_can_send_sky(start_site, sky_station, m_n);
                }else if(i>0){
                    most_close = find_most_close(start_site, can_send, start_time, skys[i-1],m_n,can_send_end);
                }else{
                    most_close = find_most_close(start_site, can_send, start_time, end_site,m_n,can_send_end);
                }
                return most_close;
            }else{
                vector<float> most_close = find_most_close(start_site, can_send, start_time, sky_station,m_n,can_send_end);
                return most_close;
            }
        }
    }
    vector<float> most_close = find_most_close(start_site, can_send, start_time, end_site,m_n,can_send_end);
    res =  most_close;


    if(can_send_end.size()!=0){
        vector<float> m_n_temp =  close_to_end(start_site, can_send_end, end_site, start_time, end_time, m_n);
        res = m_n_temp;
    }
    return res;
}

vector<float> get_next_v4(vector<vector<float>> &skys, vector<float> &start_site,vector<float> &end_site, vector<vector<int>> &can_send, float start_time, float end_time, vector<int> &m_n, vector<float> &turn_site, unordered_map<string, vector<pair<float, float>>> &record_time){
    vector<vector<int>> can_send_end;
    vector<float> res;
    for(int i=skys.size()-1;i>=0;i--){
        vector<float> sky_station = skys[i];
        if(check_inside(start_site, sky_station, end_site)){
            // cout<<compute_distance_site(start_site, sky_station)<<endl;
            // show_site(start_site);
            // show_site(sky_station);
            // int a = 
            if(compute_distance_site(start_site, sky_station) <= d){
                vector<float> most_close = check_can_send_sky(start_site, sky_station, m_n);
                // show_site(most_close);
                return most_close;
            }else{
                vector<float> most_close = find_most_close_v2(start_site, can_send, start_time, sky_station,m_n,can_send_end,record_time);
                return most_close;
            }
        }
    }
    vector<float> most_close = find_most_close_v2(start_site, can_send, start_time, end_site,m_n,can_send_end, record_time);
    res =  most_close;


    if(can_send_end.size()!=0){
        vector<float> m_n_temp =  close_to_end(start_site, can_send_end, end_site, start_time, end_time, m_n);
        res = m_n_temp;
    }
    return res;
}


float round(float number, unsigned int bits) {
    stringstream ss;
    ss << fixed << setprecision(bits) << number;
    ss >> number;
    return number;
}


string end(vector<vector<float>> &path){
    string res = "";
    for(auto site:path){
        res+="("+to_string(site[2])+","+to_string((int)site[0])+","+to_string((int)site[1])+"),";
    }
    return res.substr(0, res.size()-1);
}

bool is_intersect(float left1, float right1, float left2, float right2){
    if(right1 <= right2 && right1 >= left2 || right2 >= left1 && right2 <= right1){
        return true;
    }
    return false;
}
float check_time_is_legal(float start_time, float end_time, vector<float> &start_site, vector<float> &end_site){
    if(end_time - start_time < tf + compute_distance_site(start_site, end_site)/10000.0){
        // cout<<"***********************************************88"<<endl;
        // cout<<tf + compute_distance_site(start_site, end_site)/10000.0 - (end_time - start_time)<<endl;
        // float num = 
        // float dis = compute_distance_site(start_site, end_site)/10000.0;
        // float num = round(dis, 4) > dis ? round(dis, 4)+0.0001:round(dis, 4)+0.0001;
        return start_time + tf + compute_distance_site(start_site, end_site)/10000.0 + 0.0001;
    }
    
    return end_time+0.0001;
}

float check_time(float link_start_time,float cost_time, vector<pair<float, float>> &used_times){
    //return used time
    sort(used_times.begin(), used_times.end());
 
    for(int i=0;i<used_times.size();i++){
        pair<float,float> used_time  = used_times[i];
        if(is_intersect(link_start_time, link_start_time+cost_time, used_time.first, used_time.second)){
            if(i==used_times.size()-1 || !is_intersect(used_time.second, used_time.second+cost_time, used_times[i+1].first, used_times[i+1].second)){

                return used_time.second + cost_time;
            }
            link_start_time = used_times[i+1].second;
        }
    }
    return link_start_time + cost_time;
    // cout<<endl;
}


void get_one_path_alone_up(vector<float> &start_site_, vector<float> &end_site,vector<vector<float>>  &sky_ups, float start_time_, vector<vector<float>> &path, vector<int> &mn_,unordered_map<string, vector<pair<float, float>>> &record_time,vector<float> &ABC){
    string mode = "up";
    vector<float> turn_site = get_site(start_site_,end_site,mode, ABC);
    
    float start_time = start_time_;
    bool flag = true;
    vector<float> start_site = start_site_;
    vector<int> mn = mn_;

    while(flag){
        vector<vector<int>> can_send = find_can_send_v2(start_site, start_time, mn);
        vector<int> mn_temp = mn;
        float end_time = start_time + tf + d/10000.0;
        vector<float> start_site_temp = start_site; 

        start_site = get_next_v4(sky_ups, start_site, end_site, can_send, start_time, end_time, mn, turn_site, record_time);

        float start_time_temp = start_time;
        float cost_time= tf + compute_distance_site(start_site, start_site_temp)/10000.0;
        string query = to_string((int)mn_temp[0])+to_string((int)mn_temp[1])+to_string((int)mn[0])+to_string((int)mn[1]);
        start_time = check_time(start_time, cost_time, record_time[query]);
        start_time = check_time_is_legal(start_time_temp, start_time, start_site, start_site_temp);
        path.push_back({(float)mn[0], (float)mn[1], start_time});
        // show_site(start_site);
        if(compute_distance_site(start_site, end_site)<=D){

            flag=false;
            start_time_temp = start_time;
            float cost_time = tf + compute_distance_site(start_site, end_site)/10000.0; 
            string query = to_string((int)mn[0])+to_string((int)mn[1])+to_string((int)-100)+to_string((int)-100);
            start_time = check_time(start_time, cost_time, record_time[query]);
            start_time = check_time_is_legal(start_time_temp, start_time, start_site, end_site);
            path.push_back({-100, -100, start_time});
        }
    }
}

void get_one_path_alone_down(vector<float> &start_site_, vector<float> &end_site,vector<vector<float>>  &sky_downs, float start_time_, vector<vector<float>> &path, vector<int> &mn_,unordered_map<string, vector<pair<float, float>>> &record_time,vector<float> &ABC){
    string mode = "down";
    vector<float> turn_site = get_site(start_site_,end_site,mode, ABC);
    
    float start_time = start_time_;
    bool flag = true;
    vector<float> start_site = start_site_;
    vector<int> mn = mn_;

    while(flag){
        vector<vector<int>> can_send = find_can_send_v2(start_site, start_time, mn);
        vector<int> mn_temp = mn;
        float end_time = start_time + tf + d/10000.0;
        vector<float> start_site_temp = start_site; 

        start_site = get_next_v4(sky_downs, start_site, end_site, can_send, start_time, end_time, mn, turn_site, record_time);

        float start_time_temp = start_time;
        float cost_time= tf + compute_distance_site(start_site, start_site_temp)/10000.0;
        string query = to_string((int)mn_temp[0])+to_string((int)mn_temp[1])+to_string((int)mn[0])+to_string((int)mn[1]);
        start_time = check_time(start_time, cost_time, record_time[query]);
        start_time = check_time_is_legal(start_time_temp, start_time, start_site, start_site_temp);
        path.push_back({(float)mn[0], (float)mn[1], start_time});
        // show_site(start_site);
        if(compute_distance_site(start_site, end_site)<=D){

            flag=false;
            start_time_temp = start_time;
            float cost_time = tf + compute_distance_site(start_site, end_site)/10000.0; 
            string query = to_string((int)mn[0])+to_string((int)mn[1])+to_string((int)-100)+to_string((int)-100);
            start_time = check_time(start_time, cost_time, record_time[query]);
            start_time = check_time_is_legal(start_time_temp, start_time, start_site, end_site);
            path.push_back({-100, -100, start_time});
        }
    }
}

void get_one_path_up(vector<float> &start_site_, vector<float> &end_site,vector<vector<float>>  &sky_ups, float start_time_, vector<vector<float>> &path, vector<int> &mn_,unordered_map<string, vector<pair<float, float>>> &record_time){
    float start_time = start_time_;
    bool flag = true;
    vector<float> start_site = start_site_;
    vector<int> mn = mn_;

    while(flag){
        vector<vector<int>> can_send = find_can_send_v2(start_site, start_time, mn);
        vector<int> mn_temp = mn;
        float end_time = start_time + tf + d/10000.0;
        vector<float> start_site_temp = start_site; 

        start_site = get_next_v3(sky_ups, start_site, end_site, can_send, start_time, end_time, mn);

        float start_time_temp = start_time;
        float cost_time= tf + compute_distance_site(start_site, start_site_temp)/10000.0;
        string query = to_string((int)mn_temp[0])+to_string((int)mn_temp[1])+to_string((int)mn[0])+to_string((int)mn[1]);
        start_time = check_time(start_time, cost_time, record_time[query]);
        start_time = check_time_is_legal(start_time_temp, start_time, start_site, start_site_temp);
        path.push_back({(float)mn[0], (float)mn[1], start_time});
        // show_site(start_site);
        if(compute_distance_site(start_site, end_site)<=D){

            flag=false;
            start_time_temp = start_time;
            float cost_time = tf + compute_distance_site(start_site, end_site)/10000.0; 
            string query = to_string((int)mn[0])+to_string((int)mn[1])+to_string((int)-100)+to_string((int)-100);
            start_time = check_time(start_time, cost_time, record_time[query]);
            start_time = check_time_is_legal(start_time_temp, start_time, start_site, end_site);
            path.push_back({-100, -100, start_time});
        }
    }
}


void get_one_path_down(vector<float> &start_site_,vector<float> &end_site, vector<vector<float>>  &sky_downs, float start_time_, vector<vector<float>> &path, vector<int> &mn_,unordered_map<string, vector<pair<float, float>>> &record_time){
    float start_time = start_time_;
    bool flag = true;
    vector<float> start_site = start_site_;
    vector<int> mn = mn_;

    while(flag){
        vector<vector<int>> can_send = find_can_send_v2(start_site, start_time, mn);
        vector<int> mn_temp = mn;
        float end_time = start_time + tf + d/10000.0;
        vector<float> start_site_temp = start_site; 

        start_site = get_next_v3(sky_downs, start_site, end_site, can_send, start_time, end_time, mn);

        float start_time_temp = start_time;
        float cost_time= tf + compute_distance_site(start_site, start_site_temp)/10000.0;
        string query = to_string((int)mn_temp[0])+to_string((int)mn_temp[1])+to_string((int)mn[0])+to_string((int)mn[1]);
        start_time = check_time(start_time, cost_time, record_time[query]);
        start_time = check_time_is_legal(start_time_temp, start_time, start_site, start_site_temp);
        path.push_back({(float)mn[0], (float)mn[1], start_time});
        // show_site(start_site);
        if(compute_distance_site(start_site, end_site)<=D){

            flag=false;
            start_time_temp = start_time;
            float cost_time = tf + compute_distance_site(start_site, end_site)/10000.0; 
            string query = to_string((int)mn[0])+to_string((int)mn[1])+to_string((int)-100)+to_string((int)-100);
            start_time = check_time(start_time, cost_time, record_time[query]);
            start_time = check_time_is_legal(start_time_temp, start_time, start_site, end_site);
            path.push_back({-100, -100, start_time});
        }
    }
}

vector<string> choice_path(float start_time, vector<vector<float>> &path_up_alone, vector<vector<float>> &path_down_alone, vector<vector<float>> &path_up, vector<vector<float>> &path_down, int i_, int j_, int num_sig, unordered_map<string, vector<pair<float, float>>> &record_time){
    float time1 = path_up[path_up.size()-1][2], time2 = path_down[path_down.size()-1][2];
    float time3 = path_up_alone[path_up_alone.size()-1][2], time4 = path_down_alone[path_down_alone.size()-1][2];
    vector<vector<float>> path1 = time1>time2?path_down:path_up;
    vector<vector<float>> path2 = time3>time4?path_down_alone:path_up_alone;
    float time5 = path1[path1.size()-1][2], time6 = path2[path2.size()-1][2];
    vector<vector<float>> path = time5>time6?path2:path1;

    float time = time5>time6?time6:time5;
    string str1="",str2="",query;
    str1 = to_string(start_time)+","+to_string(i_)+","+to_string(j_)+","+to_string(time-start_time)+","+to_string(num_sig);

    if(path[0][0]>-50 && path[0][0]<50){
        str2+="("+to_string(path[0][2])+","+to_string((int)path[0][0])+","+to_string((int)path[0][1])+"),";
    }else if(path[0][0]>100){
        str2+="("+to_string(path[0][2])+","+to_string((int)(path[0][0]-100))+"),";
    }else{

    }

    for(int i=1;i<path.size();i++){
        vector<float> path_last = path[i-1];
        vector<float> path_site = path[i];
        query = to_string((int)path_last[0])+to_string((int)path_last[1])+to_string((int)path_site[0])+to_string((int)path_site[1]);
        // cout<<query<<endl;
        record_time[query].push_back({path_last[2], path_site[2]});
        if(path_site[0]>-50 && path_site[0]<50){
            str2+="("+to_string(path_site[2])+","+to_string((int)path_site[0])+","+to_string((int)path_site[1])+"),";
        }else if(path_site[0]>=100){
            str2+="("+to_string(path_site[2])+","+to_string((int)(path_site[0]-100))+"),";
        }else{
            continue;
        }

    }
    
    vector<string> res = {str1, str2.substr(0, str2.size()-1)};
    return res;
}

vector<string> plan_path_v2(float time,vector<float> &start_site,vector<float> &end_site, vector<float> &ABC, int i_, int j_){
    vector<string> res;
    float time_temp = time;
    vector<float> num_sigs = {3,3,3,1};
    vector<float> start_site_ = start_site;
    vector<vector<int>> firsts = find_first(time, D, start_site, end_site); 

    int nums = firsts.size();
    if(nums==1){
        firsts.push_back(firsts[0]);
        firsts.push_back(firsts[0]);
        firsts.push_back(firsts[0]);
    }else if(nums==2){
        firsts.push_back(firsts[0]);
        firsts.push_back(firsts[1]);
    }else if(nums==3){
        firsts.push_back(firsts[0]);
    }else if(nums==4){
        firsts[3] = firsts[0];
    }

    vector<vector<float>>  sky_ups = find_sky_station_up(ABC, start_site_, end_site);
    vector<vector<float>>  sky_downs = find_sky_station_down(ABC, start_site_, end_site);

    unordered_map<string, vector<pair<float, float>>> record_time;

   for(int i=0;i<num_sigs.size();i++){
        vector<int> first= firsts[i];
        string query = to_string(first[0])+to_string(first[1]);

        float start_time = time;

        vector<vector<float>> path_up, path_down, path_up_alone, path_down_alone;
        
        vector<float> start = {start_time*v+ first[0]*Dintra , Dinter*first[1],H};
        float s_first = compute_distance_site(start_site_, start);
        float cost_time = s_first/10000.0+tf; 
        float start_time_temp = start_time;

        start_time = check_time(start_time, cost_time, record_time[query]); //check time
        start_time = check_time_is_legal(start_time_temp, start_time, start_site_, start);

        record_time[query].push_back({start_time-cost_time, start_time});

        vector<float> start_site = {start_time*v+ first[0]*Dintra , Dinter*first[1], H};
        vector<int> mn = first;
        path_up.push_back({(float)mn[0], (float)mn[1], start_time});
        path_down.push_back({(float)mn[0], (float)mn[1], start_time});
        path_up_alone.push_back({(float)mn[0], (float)mn[1], start_time});
        path_down_alone.push_back({(float)mn[0], (float)mn[1], start_time});

        vector<int> mn_up = mn, mn_down = mn, mn_up_alone = mn, mn_down_alone = mn;

        get_one_path_up(start_site, end_site, sky_ups, start_time, path_up, mn_up, record_time);
        get_one_path_down(start_site, end_site, sky_downs, start_time, path_down, mn_down, record_time);
        get_one_path_alone_up(start_site, end_site, sky_ups, start_time, path_up_alone, mn_up_alone,record_time,ABC);
        get_one_path_alone_down(start_site, end_site, sky_downs, start_time, path_down_alone, mn_down_alone,record_time,ABC);
        vector<string> res_one_path = choice_path(time, path_up_alone, path_down_alone, path_up, path_down, i_, j_,num_sigs[i], record_time);
        res.push_back(res_one_path[0]);
        res.push_back(res_one_path[1]);
    }

    return res;
}


void done(vector<string> &one_result){
    ofstream ofs;
    ofs.open("./result.txt", ios::app);
    for(auto line:one_result){
        ofs<<line<<endl;
        cout<<line<<endl;
    }
}

int main(int argc, char const *argv[]){
    vector<float> start_time = {0, 4.7, 16.4};
    vector<vector<float>> sites = {{45.73, 45.26, 0.0}, {1200, 700, 0}, {-940, 1100, 0}};
    vector<string> res;
    for(float time:start_time){
        for(int i=0;i<sites.size();i++){
            vector<float> start_site = sites[i];
            for(int j=0;j<sites.size();j++){
                vector<float>  end_site = sites[j];
                if(i != j ){
                    vector<float> ABC = compute_ABC(start_site, end_site);
                    vector<string> one_result = plan_path_v2(time,start_site, end_site, ABC, i, j);
                    for(auto str1:one_result){
                        res.push_back(str1);
                    }
                
                }
            }

        }
    }
    done(res);
    
}
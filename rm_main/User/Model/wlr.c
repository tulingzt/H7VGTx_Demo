#include "wlr.h"
#include "chassis_task.h"
#include "leg_vmc.h"
#include "wheel_leg_model.h"
#include "prot_imu.h"
#include "drv_dji_motor.h"
#include "pid.h"
#include "kalman_filter.h"
#include "math_lib.h"
#include "math_matrix.h"

#define WLR_SIGN(x) ((x) > 0? (1): (-1))

#define CHASSIS_PERIOD_DU 2

const float LegLengthParam[5] = {0.150f, 0.270f, 0.270f, 0.150f, 0.150f};
float mb = 4.4f, ml = 2.09f, mw = 0.715f;//机体质量 腿部质量 轮子质量 14.5
const float BodyWidth = 0.51f;//两轮间距
//const float WheelRadius = 0.10f;//轮子半径
const float WheelRadius = 0.065f;//轮子半径
const float LegLengthMax = 0.35f, LegLengthMin = 0.11f;

const float LegLengthJump1 = 0.12f;//压腿
const float LegLengthJump2 = 0.35f;//蹬腿
const float LegLengthJump3 = 0.18f;//收腿
const float LegLengthJump4 = 0.15f;//落地

const float LegLengthHightFly = 0.30f;//长腿腿长腾空 0.28
const float LegLengthFly = 0.25f;//正常腿长腾空
const float LegLengthHigh = 0.20f;//长腿 0.23
const float LegLengthNormal = 0.15f;//正常

//const float LegLengthHightFly = 0.28f;//长腿腿长腾空 0.28
//const float LegLengthFly = 0.20f;//正常腿长腾空
//const float LegLengthHigh = 0.23f;//长腿 0.23
//const float LegLengthNormal = 0.12f;//正常

float x3_balance_zero = 0.00f, x5_balance_zero = 0.00f;//腿摆角角度偏置 机体俯仰角度偏置
float x3_fight_zero = 0.06f;
//								位移  速度	角度	角速度  角度	角速度
float K_Array_Wheel[2][6] =		{{60, 30, 80, 8, 300, 10}, 
                                { -0, -0.7, -8, -1, 3, 2}};
float K_Array_Leg[2][6] =		{{6.08535, 12.6636, 49.1418, 7.57203, 54.8073, 12.7387}, {-0.793527, -1.75976, -13.1544, -1.78789, 4.8796, 1.66868}};
float K_Array_Fly[2][6] =		{{0, 0, 80, 10, 0, 0}, { 0, 0, 0, 0, 0, 0}};
float K_Array_Prone[2][6] =     {{0, 0, 0, 0, 0, 0}, {0, 20, 0, 0, 0, 2}};

float K_Array_test[2][6];
float K_Array_List[12][4] = 
//{{12.1359,15.7251,-134.375,159.136},{24.9973,4.83081,-172.951,216.835},{36.757,497.27,-1734.71,1671.41},{11.3919,38.5109,-142.095,128.464},{2.83296,166.971,228.572,-623.115},{-1.87233,27.389,38.6001,-101.161},{-0.837065,-0.595231,-15.3135,25.0277},{-1.79608,0.294299,-32.6456,49.9669},{-5.80305,-73.8945,-20.7886,112.117},{-0.78318,-2.60947,-31.1545,37.5597},{15.2539,-1.90842,-76.2986,95.478},{2.92341,-2.30254,-4.92724,7.0132}}; 
    //低可以 高不行
//{{12.5193,14.6442,-124.321,146.563},{36.3924,25.8018,-294.012,351.744},{38.9382,584.8,-1843.53,1700.3},{9.73159,56.8841,-135.027,95.889},{-5.12791,180.448,215.719,-615.286},{-4.03071,31.2338,36.1287,-101.436},{-0.699727,-0.799938,-15.1699,24.9683},{-2.1674,-1.09018,-46.6402,73.3611},{-6.17683,-71.0743,-73.3065,168.045},{-0.841997,-1.12008,-44.9368,48.7584},{15.6952,-2.62929,-65.4861,80.7837},{3.03518,-2.4289,-2.61137,3.7334}};

//{{13.4634,-19.7832,-54.5177,113.387},{24.4316,-51.3459,-32.3451,127.275},{48.7587,185.387,-1117.15,1389.19},{10.8177,7.18266,-104.476,142.102},{8.8817,529.301,-1060.54,729.937},{0.475421,112.637,-221.565,149.593},{-0.747635,-7.82155,13.8973,-7.59682},{-1.53314,-11.2506,17.1682,-7.14969},{-4.36387,-100.154,149.184,-87.7787},{-0.578323,-8.85429,-1.62092,6.60359},{19.7112,-36.7279,-31.8389,101.513},{4.28528,-7.70948,-5.93745,19.9717}};
{{13.4664,-20.2539,-54.009,113.604},{23.1388,-51.262,-22.0911,111.965},{48.162,174.661,-1084.4,1358.26},{10.8266,4.86413,-99.6066,138.581},{9.6221,529.676,-1068.49,740.669},{0.688052,112.553,-223.065,151.83},{-0.751659,-7.88392,14.1798,-7.917},{-1.47091,-10.4386,16.0673,-6.81278},{-4.30545,-99.3549,152.465,-92.1674},{-0.571589,-8.8386,-0.292307,5.31904},{19.6674,-37.2007,-31.5412,102.148},{4.27126,-7.81527,-5.88522,20.1278}};
    float K_Array_List_lq[12][4][4] = 
//{{{6.13556,-0.00812017,0.000047346,-7.85689e-8},{11.4091,0.0147766,3.73834e-6,0},{-123.106,-0.0343755,0,0},{202.509,0,0,0}},
//{{16.8765,-0.0215183,0.000126153,-2.13668e-7},{22.3487,0.0414594,-1.65244e-6,0},{-290.877,-0.0918103,0,0},{486.196,0,0,0}},
//{{12.0713,-0.0357215,0.000187603,-1.5406e-7},{468.34,-0.0211163,0.000421664,0},{-1981.22,-0.11515,0,0},{2561.85,0,0,0}},
//{{3.35549,-0.000952253,6.19623e-6,-1.70836e-8},{26.1622,0.00161176,7.45819e-6,0},{-81.3821,-0.00749412,0,0},{77.9294,0,0,0}},
//{{1.35547,0.0837986,-0.000494484,8.91191e-7},{116.216,-0.23393,0.00035451,0},{282.336,0.384397,0,0},{-923.119,0,0,0}},
//{{-0.387073,0.0199651,-0.000118316,2.16491e-7},{27.0777,-0.0564747,0.0000873809,0},{68.6347,0.0922291,0,0},{-218.658,0,0,0}},
//{{-0.564074,-0.00239842,0.0000139335,-2.38582e-8},{-0.691714,0.00625592,-7.94146e-6,0},{-17.3956,-0.010863,0,0},{38.8056,0,0,0}},
//{{-1.69707,-0.00660206,0.0000382862,-6.53947e-8},{-0.046062,0.0160347,-0.0000152172,0},{-54.856,-0.0297814,0,0},{114.731,0,0,0}},
//{{-3.57855,-0.0281303,0.000172425,-3.33176e-7},{-92.4118,0.0818674,-0.000149969,0},{-14.7022,-0.125051,0,0},{195.299,0,0,0}},
//{{-0.80101,-0.00510864,0.0000304449,-5.38814e-8},{3.55887,0.00312333,0.0000348456,0},{-66.519,-0.0198739,0,0},{99.6118,0,0,0}},
//{{13.8316,-0.00881745,0.0000534308,-9.54403e-8},{-6.99012,0.0213533,-0.0000292566,0},{-89.2701,-0.036652,0,0},{161.716,0,0,0}},
//{{3.32847,-0.00181106,0.0000111953,-2.13373e-8},{-2.22601,0.00514635,-9.72938e-6,0},{-15.8665,-0.00787378,0,0},{29.0107,0,0,0}}};

{{{11.5204,-0.0182899,0.000106213,-1.66495e-7},{31.3417,0.030098,0.0000151576,0},{-296.919,-0.0724147,0,0},{482.446,0,0,0}},
{{24.0908,-0.035299,0.000207037,-3.38145e-7},{23.0936,0.0664406,-0.0000130899,0},{-411.894,-0.141962,0,0},{711.923,0,0,0}},
{{26.2754,-0.0811866,0.000444932,-4.73673e-7},{672.615,0.0121388,0.000618139,0},{-3174.89,-0.26074,0,0},{4291.53,0,0,0}},
{{9.55257,-0.0120772,0.0000714574,-1.18146e-7},{43.1285,0.0186104,0.0000169431,0},{-247.237,-0.0483437,0,0},{342.838,0,0,0}},
{{15.6127,0.188081,-0.00109772,1.80116e-6},{170.75,-0.504943,0.000846888,0},{911.139,0.792457,0,0},{-2429.15,0,0,0}},
{{1.76811,0.0400676,-0.000234735,3.91483e-7},{36.8213,-0.10831,0.000181866,0},{194.31,0.170112,0,0},{-514.487,0,0,0}},
{{-0.96691,-0.00351856,0.0000202789,-3.1637e-8},{0.27558,0.00897766,-0.0000136596,0},{-30.1037,-0.0146195,0,0},{63.1165,0,0,0}},
{{-2.13968,-0.00658082,0.0000377782,-5.90553e-8},{3.86548,0.0150992,-0.0000151663,0},{-68.5344,-0.0275499,0,0},{132.807,0,0,0}},
{{-4.35166,-0.0298867,0.000181951,-3.36782e-7},{-92.7861,0.087823,-0.000176006,0},{-20.3266,-0.127832,0,0},{221.012,0,0,0}},
{{-0.963148,-0.00576085,0.0000341389,-5.87509e-8},{2.15582,0.0057035,0.000026582,0},{-64.6508,-0.0223969,0,0},{104.1,0,0,0}},
{{18.569,-0.0178054,0.000105747,-1.74543e-7},{7.36847,0.0347009,-0.000019059,0},{-230.788,-0.0697888,0,0},{392.631,0,0,0}},
{{4.15806,-0.00335669,0.0000202361,-3.53769e-8},{-0.0427919,0.00741395,-7.64787e-6,0},{-39.3092,-0.0136382,0,0},{67.4881,0,0,0}}};

wlr_t wlr;
lqr_t lqr[2];

kalman_filter_t kal_fn[2], kal_v[2];

pid_t pid_leg_length[2];
pid_t pid_leg_length_fast[2];
pid_t pid_q0, pid_roll, pid_wx, pid_yaw, pid_wz;

static float wlr_fn_calc(float az, float Fy_fdb, float T0_fdb, float L0[3], float theta[3])
{
    float Fwy = Fy_fdb * cosf(theta[0]) + T0_fdb * sinf(theta[0]) / L0[0];//轮子受到腿部机构竖直方向的作用力
    float yw_ddot = az
                    - L0[2] * cosf(theta[0])
                    + 2 * L0[1] * theta[1] * sinf(theta[0])
                    + L0[0] * theta[2] * sinf(theta[0])
                    + L0[0] * powf(theta[1], 2) * cosf(theta[0]);//轮子竖直方向的加速度
    return Fwy + mw * GRAVITY + mw * yw_ddot;
}

static void k_array_fit1(float K[2][6], float high_fdb)
{
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 6; j++) {
            K[i][j] = 0;
            for(int h = 0; h < 4; h++) {
                K[i][j] += K_Array_List[i * 6 + j][h] * powf(high_fdb, h);
            }
        }
    }
}

static void k_array_fit2(float K[2][6], float high_fdb, float q0_fdb)
{
    float temp;
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 6; j++) {
            temp = 0;
            for(int x = 0; x <= 3; x++)
                for(int y = 0; y <= 3; y++)
                    temp += (K_Array_List_lq[i * 6 + j][x][y] * powf(high_fdb, x) * powf(q0_fdb, y));
            K[i][j] = temp;
        }
}

void wlr_init(void)
{
	wlr.high_set = LegLengthNormal;
	wlr.roll_set = 0;
	wlr.q0_set = PI / 2 + x3_balance_zero;
	
	twm_init(&twm, BodyWidth, WheelRadius);
	tlm_init(&tlm, LegLengthMax, LegLengthMin, BodyWidth);
	for(int i = 0; i < 2; i++)
	{
		//腿部长度初始化
		vmc_init(&vmc[i], LegLengthParam);
		//卡尔曼滤波器初始化
        kalman_filter_init(&kal_fn[i], 1, 0, 1);
        kal_fn[i].A_data[0] = 1;
        kal_fn[i].H_data[0] = 1;
        kal_fn[i].Q_data[0] = 1;
        kal_fn[i].R_data[0] = 100;
        kalman_filter_init(&kal_v[i], 2, 0, 2);
        kal_v[i].A_data[0] = 1; kal_v[i].A_data[1] = 0.002f; kal_v[i].A_data[3] = 1;
        kal_v[i].H_data[0] = 1; kal_v[i].H_data[3] = 1;
        kal_v[i].Q_data[0] = 0.01f; kal_v[i].Q_data[3] = 0.01f;
        kal_v[i].R_data[0] = 1.0f; kal_v[i].R_data[3] = 10.0f;
		//PID参数初始化
//		pid_init(&pid_leg_length[i], NONE, 800, 2.5f, 30000, 25, 40);//i 2.5f
        pid_init(&pid_leg_length[i], NONE, 500, 0.0f, 10000, 25, 40);//i 2.5f
		pid_init(&pid_leg_length_fast[i], NONE, 1000, 0, 10000, 0, 50);
	}
	//卡尔曼滤波器初始化

	//PID参数初始化
    pid_init(&pid_yaw, NONE, -5, 0, 0, 0, 10);
	pid_init(&pid_wz, NONE, 2.0f, 0, 7.0f, 0, 3.0f);
	pid_init(&pid_q0, NONE, 60, 0, 100, 0, 10);//与LQR的虚拟腿摆角控制拮抗 60 0 120
	pid_init(&pid_roll, NONE, 500, 0, 3000, 0, 30);//与VMC的腿长控制协同  1000 0 3500
}

void wlr_protest(void)
{
	pid_leg_length[0].i_out = 0;
	pid_leg_length[1].i_out = 0;
}

//轮子：位移、速度   摆角：角度、角速度   机体俯仰：角度、角速度
void wlr_control(void)
{
	//------------------------反馈数据更新------------------------//
    //限制
    float v_fdb = (-wlr.side[0].wy * WheelRadius-wlr.side[1].wy * WheelRadius)/2.0f;
    if (fabs(v_fdb) > fabs(wlr.v_set))//加强超速控制
        wlr.v_set = data_fusion(wlr.v_set, 0, fabs(v_fdb - wlr.v_set));
    if (chassis.mode == CHASSIS_MODE_REMOTER_ROTATE ||
        chassis.mode == CHASSIS_MODE_KEYBOARD_ROTATE) {
        wlr.wz_set = 9;                                                 //小陀螺时固定转速 减小前进
        wlr.v_set = -0.5f * wlr.v_set;
    } else {
        wlr.wz_set = pid_calc(&pid_yaw, wlr.yaw_fdb + wlr.yaw_err, wlr.yaw_fdb);
        wlr.v_set = data_fusion(wlr.v_set, 0, fabs(wlr.wz_fdb)/6.0f);   //限制前进速度 旋转过快
//        wlr.v_set = data_fusion(wlr.v_set, 0, fabs(wlr.pit_fdb)/0.3f);   //限制前进速度 底盘pit太偏
    }
    wlr.wz_set = data_fusion(wlr.wz_set, 0, fabs(v_fdb/3.0f));          //限制旋转速度 前进过快
    wlr.wz_set = data_fusion(wlr.wz_set, 0, fabs(wlr.roll_fdb/0.2f));   //限制旋转速度 roll太偏
	//更新两轮模型
	twm_feedback_calc(&twm, wlr.side[0].wy, wlr.side[1].wy, wlr.wz_fdb);//输入左右轮子转速
	twm_reference_calc(&twm, wlr.v_set, wlr.wz_set);//计算两侧轮腿模型的设定速度
	//两侧轮腿分别更新数据
	for(int i = 0; i < 2; i++) {
		//更新腿部VMC模型
		vmc_forward_solution(&vmc[i], wlr.side[i].q1, wlr.side[i].q4, wlr.side[i].w1, \
									  wlr.side[i].w4, wlr.side[i].t1, wlr.side[i].t4);
		//LQR输入反馈值
		lqr[i].last_x2 = lqr[i].X_fdb[1];
        
        kal_v[i].measured_vector[0] = -wlr.side[i].wy * WheelRadius;
        kal_v[i].measured_vector[1] = -chassis_imu.ax;
        kalman_filter_update(&kal_v[i]);
        wlr.side[i].v_fdb = -wlr.side[i].wy * WheelRadius;
        wlr.side[i].a_fdb = -chassis_imu.ax;
        wlr.side[i].v_kal = kal_v[i].filter_vector[0];
        wlr.side[i].a_kal = kal_v[i].filter_vector[1];

        lqr[i].X_fdb[1] = kal_v[i].filter_vector[0];
//		lqr[i].X_fdb[1] = -wlr.side[i].wy * WheelRadius;
        lqr[i].X_fdb[0] = -2 * WLR_SIGN(i) * driver_motor[i].position * WheelRadius;
		lqr[i].X_fdb[4] = x5_balance_zero + wlr.pit_fdb;
		lqr[i].X_fdb[5] = wlr.wy_fdb;
        if (chassis.mode == CHASSIS_MODE_KEYBOARD_FIGHT)
            lqr[i].X_fdb[2] = x3_fight_zero + (PI / 2 - lqr[i].X_fdb[4] - vmc[i].q_fdb[0]);
        else
            lqr[i].X_fdb[2] = x3_balance_zero + (PI / 2 - lqr[i].X_fdb[4] - vmc[i].q_fdb[0]);
		lqr[i].X_fdb[3] = lqr[i].X_fdb[5] - vmc[i].V_fdb.e.w0_fdb;
		lqr[i].dot_x4 = (lqr[i].X_fdb[3] - lqr[i].last_x4) / (CHASSIS_PERIOD_DU * 0.001f); //腿倾角加速度(状态变量x4的dot)计算
		lqr[i].last_x4 = lqr[i].X_fdb[3];
        if(ABS(wlr.v_set) > 1e-3f || ABS(wlr.wz_set) > 0.1f || ABS(vmc[0].q_fdb[0] - vmc[1].q_fdb[0]) > 0.05f) {//有输入速度 或 两腿有差角时 将位移反馈置0  不发挥作用
            lqr[i].X_ref[0] = lqr[i].X_fdb[0];
        }
        data_limit(&lqr[i].X_ref[0], lqr[i].X_fdb[0]-0.5f, lqr[i].X_fdb[0]+0.5f);
		//支持力解算
		float L0_array[3] = {vmc[i].L_fdb, vmc[i].V_fdb.e.vy0_fdb, vmc[i].Acc_fdb.L0_ddot};
		float theta_array[3] = {lqr[i].X_fdb[2], lqr[i].X_fdb[3], lqr[i].dot_x4};
		wlr.side[i].Fn_fdb = wlr_fn_calc(wlr.az_fdb, vmc[i].F_fdb.e.Fy_fdb, vmc[i].F_fdb.e.T0_fdb, L0_array, theta_array);
        kal_fn[i].measured_vector[0] = wlr.side[i].Fn_fdb;
        kalman_filter_update(&kal_fn[i]);
        wlr.side[i].Fn_kal = kal_fn[i].filter_vector[0];
		//离地检测
        if(wlr.side[i].Fn_kal < 20.0f)
            wlr.side[i].fly_cnt++;
        else if(wlr.side[i].fly_cnt > 0)
            wlr.side[i].fly_cnt-=2;
        if(wlr.side[i].fly_cnt > 30) {
            wlr.side[i].fly_cnt = 30;
            wlr.side[i].fly_flag = 1;
        } else if(wlr.side[i].fly_cnt == 0)
            wlr.side[i].fly_flag = 0;
	}
	//高度选择 跳跃状态改变
	if (wlr.jump_flag == 1) {//跳跃起跳状态 先压腿
		wlr.high_set = LegLengthJump1;
		if(fabs(vmc[0].L_fdb - LegLengthJump1) < 0.02f && fabs(vmc[1].L_fdb - LegLengthJump1) < 0.02f)
			wlr.jump_flag = 2;
	} else if (wlr.jump_flag == 2) {//起跳 弹腿
		wlr.high_set = LegLengthJump2;
		if(fabs(vmc[0].L_fdb - LegLengthJump2) < 0.02f && fabs(vmc[1].L_fdb - LegLengthJump2) < 0.02f)
			wlr.jump_flag = 3;
	} else if (wlr.jump_flag == 3) {//收腿
		wlr.high_set = LegLengthJump3;
//		if (fabs(vmc[0].L_fdb - LegLengthJump3) < 0.02f && fabs(vmc[1].L_fdb - LegLengthJump3) < 0.02f && !wlr.side[0].fly_flag && !wlr.side[1].fly_flag)
        if (fabs(vmc[0].L_fdb - LegLengthJump3) < 0.02f && fabs(vmc[1].L_fdb - LegLengthJump3) < 0.02f)
			wlr.jump_flag = 4;
	} else if (wlr.jump_flag == 4) {//落地
		wlr.high_set = LegLengthJump4;
		wlr.jump_cnt++;
		if (wlr.jump_cnt > 200) {
			wlr.jump_flag = 0;
			wlr.jump_cnt = 0;
		}
	} else if (wlr.high_flag) { //长腿长
        if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {//腾空
            wlr.high_set = LegLengthHightFly;
        } else {
            wlr.high_set = LegLengthHigh;
        }
    } else { //正常腿长
        if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {//腾空
            wlr.high_set = LegLengthFly;
        } else {
            wlr.high_set = LegLengthNormal;
        }
    }
    //限制 旋转压腿长用于解决侧翻 腿摆压腿长用于解决膝关节可能卡地面
//    wlr.high_set = data_fusion(wlr.high_set, LegLengthMin, (fabs(lqr[0].X_fdb[2])+fabs(lqr[1].X_fdb[2]))/2.0f/2.0f);
//     wlr.high_set = data_fusion(wlr.high_set, wlr.high_set-0.1f, fabs(wlr.pit_fdb)/0.3f);
//    wlr.high_set = wlr.high_set/arm_cos_f32((fabs(lqr[0].X_fdb[2])+fabs(lqr[1].X_fdb[2]))/2.0f);
	//更新两腿模型
	tlm_gnd_roll_calc(&tlm, -wlr.roll_fdb, vmc[0].L_fdb, vmc[1].L_fdb);//计算地形倾角
	if (wlr.jump_flag != 0 || (wlr.side[0].fly_flag && wlr.side[1].fly_flag))
		tlm.l_ref[0] = tlm.l_ref[1] = wlr.high_set;
	else
        tlm_leg_length_cacl(&tlm, wlr.high_set, 0);//计算腿长设定值
	//------------------------状态选择------------------------//
	//根据当前状态选择合适的控制矩阵
    if (wlr.ctrl_mode == 2) {//力控
        if (wlr.prone_flag) {
            aMartix_Cover(lqr[0].K, (float*)K_Array_Prone, 2, 6);
            aMartix_Cover(lqr[1].K, (float*)K_Array_Prone, 2, 6);
        } else if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {//腾空
            aMartix_Cover(lqr[0].K, (float*)K_Array_Fly, 2, 6);
            aMartix_Cover(lqr[1].K, (float*)K_Array_Fly, 2, 6);
        } else {
            k_array_fit1(K_Array_test, vmc[0].L_fdb);
            aMartix_Cover(lqr[0].K, (float*)K_Array_test, 2, 6);
            k_array_fit1(K_Array_test, vmc[1].L_fdb);
            aMartix_Cover(lqr[1].K, (float*)K_Array_test, 2, 6);
//            k_array_fit2(K_Array_test, vmc[0].L_fdb, vmc[0].q_fdb[0]/PI*180);
//            aMartix_Cover(lqr[0].K, (float*)K_Array_test, 2, 6);
//            k_array_fit2(K_Array_test, vmc[1].L_fdb, vmc[1].q_fdb[0]/PI*180);
//            aMartix_Cover(lqr[1].K, (float*)K_Array_test, 2, 6);
        }
    } else if (wlr.ctrl_mode == 1) {//位控
        if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {//腾空
            aMartix_Cover(lqr[0].K, (float*)K_Array_Fly, 2, 6);
            aMartix_Cover(lqr[1].K, (float*)K_Array_Fly, 2, 6);
        } else {
            aMartix_Cover(lqr[0].K, (float*)K_Array_Wheel, 2, 6);
            aMartix_Cover(lqr[1].K, (float*)K_Array_Wheel, 2, 6);
        }
    }
	//------------------------控制数据更新------------------------//
	//全身运动控制
    if (fabs(wlr.wz_set) > 1) {
        wlr.q0_offs   = 0;                                                  //小陀螺时不用同步
    } else {
        wlr.q0_offs = pid_calc(&pid_q0, vmc[0].q_fdb[0], vmc[1].q_fdb[0]);//双腿摆角同步控制
    }
	wlr.roll_offs = pid_calc(&pid_roll, wlr.roll_set, wlr.roll_fdb);
    if ((wlr.side[0].fly_flag && wlr.side[1].fly_flag) || fabs(wlr.pit_fdb)>0.15f) {//腾空或者未平衡时不旋转
        wlr.wz_offs   = 0;
    } else {
        wlr.wz_offs   = pid_calc(&pid_wz, wlr.wz_fdb, wlr.wz_set);//Yaw控制
    }
	//两侧轮腿分别独立控制
	for (int i = 0; i < 2; i++) {
		//LQR输入控制值 计算得出轮子与腿的力矩
        //速度融合 在快要撞限位时减小速度控制
//        if (wlr.side[i].q1 > 3.14f)
//            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (wlr.side[i].q1 - 3.14f)/0.5f);
//        else if (wlr.side[i].q4 < -0.0f)
//            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (-0.0f - wlr.side[i].q4)/0.5f);
//        else
//            lqr[i].X_ref[1] = twm.v_ref[i];
        
//        if(wlr.side[i].q1 > 3.3f)
//            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (wlr.side[i].q1 - 3.3f)/0.4f);
//        else if(wlr.side[i].q4 < -0.2f)
//            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (-0.2f - wlr.side[i].q4)/0.4f);
//        else
//            lqr[i].X_ref[1] = twm.v_ref[i];
        if(wlr.side[i].q1 > 3.6f)
            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (wlr.side[i].q1 - 3.6f)/0.1f);
        else if(wlr.side[i].q4 < -0.5f)
            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (-0.5f - wlr.side[i].q4)/0.1f);
        else
            lqr[i].X_ref[1] = twm.v_ref[i];
        
		aMartix_Add(1, lqr[i].X_ref, -1, lqr[i].X_fdb, lqr[i].X_diff, 6, 1);
		aMartix_Mul(lqr[i].K, lqr[i].X_diff, lqr[i].U_ref, 2, 6, 1);
		//腿部虚拟力控制
		if (wlr.jump_flag == 2 || wlr.jump_flag == 4)						//跳跃蹬腿阶段 响应要大
			wlr.side[i].Fy = pid_calc(&pid_leg_length_fast[i], tlm.l_ref[i], vmc[i].L_fdb)\
								+ 40.0f + WLR_SIGN(i) * wlr.roll_offs;
		else if (wlr.jump_flag == 3)                                        //跳跃收腿阶段 响应要大
			wlr.side[i].Fy = pid_calc(&pid_leg_length_fast[i], tlm.l_ref[i], vmc[i].L_fdb);
		else if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {            //浮空收腿 响应不用那么大
			wlr.side[i].Fy = pid_calc(&pid_leg_length_fast[i], tlm.l_ref[i], vmc[i].L_fdb);
		} else																//常态 跳跃压腿阶段 跳跃落地阶段
			wlr.side[i].Fy = pid_calc(&pid_leg_length[i], tlm.l_ref[i], vmc[i].L_fdb)\
                                 + 38.0f + WLR_SIGN(i) * wlr.roll_offs;
		wlr.side[i].T0 = -lqr[i].U_ref[0] / 2 + WLR_SIGN(i) * wlr.q0_offs;	//两条腿时，LQR输出力矩需要除以2
		vmc_inverse_solution(&vmc[i], wlr.high_set, wlr.q0_set, wlr.side[i].T0, wlr.side[i].Fy);
	}
	//------------------------控制数据输出------------------------//
	for (int i = 0; i < 2; i++) {
		wlr.side[i].T1 =  vmc[i].T_ref.e.T1_ref;
		wlr.side[i].T4 =  vmc[i].T_ref.e.T4_ref;
        wlr.side[i].Tw = -lqr[i].U_ref[1] + WLR_SIGN(i) * wlr.wz_offs;
		wlr.side[i].P1 =  vmc[i].q_ref[1];
		wlr.side[i].P4 =  vmc[i].q_ref[4];
	}
}

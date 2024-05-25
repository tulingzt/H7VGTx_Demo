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
const float WheelRadius = 0.065f;//轮子半径
const float LegLengthMax = 0.35f, LegLengthMin = 0.11f;

const float LegLengthJump1 = 0.12f;//压腿
const float LegLengthJump2 = 0.35f;//蹬腿
const float LegLengthJump3 = 0.18f;//收腿
const float LegLengthJump4 = 0.15f;//落地

const float LegLengthHightFly = 0.30f;//长腿腿长腾空 0.28
const float LegLengthFly = 0.25f;//正常腿长腾空
const float LegLengthHigh2 = 0.30f;//超长腿
const float LegLengthHigh = 0.20f;//长腿 0.23
const float LegLengthNormal = 0.15f;//正常

float x3_balance_zero = 0.00f, x5_balance_zero = 0.00f;//腿摆角角度偏置 机体俯仰角度偏置
float x3_fight_zero = -0.01f;
//								位移  速度	角度	角速度  角度	角速度
float K_Array_Wheel[2][6] =		{{60, 30, 80, 8, 300, 10}, 
                                { -0, -0.7, -8, -1, 3, 2}};
float K_Array_Leg[2][6] =		{{6.08535, 12.6636, 49.1418, 7.57203, 54.8073, 12.7387}, {-0.793527, -1.75976, -13.1544, -1.78789, 4.8796, 1.66868}};
float K_Array_Fly[2][6] =		{{0, 0, 80, 10, 0, 0}, { 0, 0, 0, 0, 0, 0}};
float K_Array_Prone[2][6] =     {{0, 0, 0, 0, 0, 0}, {0, 20, 0, 0, 0, 2}};

float K_Array_test[2][6];
float K_Array_List[12][4] = 
//{{13.4664,-20.2539,-54.009,113.604},{23.1388,-51.262,-22.0911,111.965},{48.162,174.661,-1084.4,1358.26},{10.8266,4.86413,-99.6066,138.581},{9.6221,529.676,-1068.49,740.669},{0.688052,112.553,-223.065,151.83},{-0.751659,-7.88392,14.1798,-7.917},{-1.47091,-10.4386,16.0673,-6.81278},{-4.30545,-99.3549,152.465,-92.1674},{-0.571589,-8.8386,-0.292307,5.31904},{19.6674,-37.2007,-31.5412,102.148},{4.27126,-7.81527,-5.88522,20.1278}};
    //加强定点和速度 150 150
{{16.5023,-23.0667,-68.527,139.098},{28.4325,-56.4708,-45.4545,153.884},{50.8603,206.621,-1188.61,1460.},{10.8696,12.4379,-116.584,151.693},{6.66386,533.585,-1057.27,719.658},{-0.149511,114.118,-221.847,148.046},{-0.898044,-9.47829,16.4422,-8.60269},{-1.76035,-13.0749,19.0238,-7.05895},{-4.53472,-102.181,142.023,-77.9829},{-0.596324,-8.95685,-4.5413,9.44794},{19.8547,-35.9135,-32.499,100.511},{4.32843,-7.53239,-6.02343,19.6762}};
    
float K_Array_List_lq[12][4][4] = 
//{{{11.5204,-0.0182899,0.000106213,-1.66495e-7},{31.3417,0.030098,0.0000151576,0},{-296.919,-0.0724147,0,0},{482.446,0,0,0}},
//{{24.0908,-0.035299,0.000207037,-3.38145e-7},{23.0936,0.0664406,-0.0000130899,0},{-411.894,-0.141962,0,0},{711.923,0,0,0}},
//{{26.2754,-0.0811866,0.000444932,-4.73673e-7},{672.615,0.0121388,0.000618139,0},{-3174.89,-0.26074,0,0},{4291.53,0,0,0}},
//{{9.55257,-0.0120772,0.0000714574,-1.18146e-7},{43.1285,0.0186104,0.0000169431,0},{-247.237,-0.0483437,0,0},{342.838,0,0,0}},
//{{15.6127,0.188081,-0.00109772,1.80116e-6},{170.75,-0.504943,0.000846888,0},{911.139,0.792457,0,0},{-2429.15,0,0,0}},
//{{1.76811,0.0400676,-0.000234735,3.91483e-7},{36.8213,-0.10831,0.000181866,0},{194.31,0.170112,0,0},{-514.487,0,0,0}},
//{{-0.96691,-0.00351856,0.0000202789,-3.1637e-8},{0.27558,0.00897766,-0.0000136596,0},{-30.1037,-0.0146195,0,0},{63.1165,0,0,0}},
//{{-2.13968,-0.00658082,0.0000377782,-5.90553e-8},{3.86548,0.0150992,-0.0000151663,0},{-68.5344,-0.0275499,0,0},{132.807,0,0,0}},
//{{-4.35166,-0.0298867,0.000181951,-3.36782e-7},{-92.7861,0.087823,-0.000176006,0},{-20.3266,-0.127832,0,0},{221.012,0,0,0}},
//{{-0.963148,-0.00576085,0.0000341389,-5.87509e-8},{2.15582,0.0057035,0.000026582,0},{-64.6508,-0.0223969,0,0},{104.1,0,0,0}},
//{{18.569,-0.0178054,0.000105747,-1.74543e-7},{7.36847,0.0347009,-0.000019059,0},{-230.788,-0.0697888,0,0},{392.631,0,0,0}},
//{{4.15806,-0.00335669,0.0000202361,-3.53769e-8},{-0.0427919,0.00741395,-7.64787e-6,0},{-39.3092,-0.0136382,0,0},{67.4881,0,0,0}}};

{{{14.1616,-0.0218561,0.000127016,-2.00385e-7},{38.2525,0.0352282,0.0000232427,0},{-358.311,-0.0869167,0,0},{579.604,0,0,0}},
{{27.8953,-0.0399102,0.000234334,-3.85462e-7},{27.8686,0.0739839,-6.58209e-6,0},{-471.999,-0.161251,0,0},{809.653,0,0,0}},
{{27.7401,-0.0814256,0.000446498,-4.81686e-7},{695.715,0.00456115,0.000668559,0},{-3240.36,-0.263723,0,0},{4351.34,0,0,0}},
{{9.63532,-0.0119904,0.0000711665,-1.19563e-7},{47.773,0.0172065,0.0000249021,0},{-253.325,-0.0483994,0,0},{341.8,0,0,0}},
{{13.725,0.189194,-0.00110557,1.82442e-6},{172.309,-0.506322,0.000838248,0},{917.462,0.799123,0,0},{-2435.82,0,0,0}},
{{1.21353,0.0404646,-0.000237418,3.98642e-7},{37.6838,-0.109153,0.000181208,0},{194.693,0.172305,0,0},{-515.166,0,0,0}},
{{-1.17392,-0.00431216,0.0000248655,-3.89244e-8},{0.536107,0.0109156,-0.0000161698,0},{-37.7648,-0.0179499,0,0},{78.3611,0,0,0}},
{{-2.46334,-0.00766015,0.0000440023,-6.899e-8},{4.80664,0.0173686,-0.0000164408,0},{-81.2759,-0.0320803,0,0},{156.149,0,0,0}},
{{-4.61797,-0.0313264,0.000190526,-3.51661e-7},{-90.8247,0.0889503,-0.000167776,0},{-43.9587,-0.133422,0,0},{254.861,0,0,0}},
{{-1.00671,-0.00600785,0.0000356169,-6.13615e-8},{2.70994,0.00563627,0.000029413,0},{-70.4931,-0.0233141,0,0},{111.428,0,0,0}},
{{18.7556,-0.0171081,0.000101849,-1.69631e-7},{6.1897,0.0331916,-0.0000170689,0},{-220.557,-0.067267,0,0},{373.99,0,0,0}},
{{4.21233,-0.00318446,0.0000192734,-3.41253e-8},{-0.379286,0.00709104,-7.47767e-6,0},{-36.6096,-0.0130026,0,0},{62.6967,0,0,0}}};

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
//		pid_init(&pid_leg_length[i], NONE, 800, 0.0f, 30000, 25, 40);//i 2.5f
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
        wlr.v_set = 0.75f * wlr.v_set;
    } else {
        wlr.wz_set = pid_calc(&pid_yaw, wlr.yaw_fdb + wlr.yaw_err, wlr.yaw_fdb);
        wlr.v_set = data_fusion(wlr.v_set, 0, fabs(wlr.wz_fdb)/6.0f);   //限制前进速度 旋转过快
        wlr.v_set = data_fusion(wlr.v_set, 0, fabs(wlr.pit_fdb)/0.3f);   //限制前进速度 底盘pit太偏
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
        if(ABS(wlr.v_set) > 1e-3f || ABS(wlr.wz_set) > 0.1f || ABS(vmc[0].q_fdb[0] - vmc[1].q_fdb[0]) > 0.07f) {//有输入速度 或 两腿有差角时 将位移反馈置0  不发挥作用
            lqr[i].X_ref[0] = lqr[i].X_fdb[0];
        }
        data_limit(&lqr[i].X_ref[0], lqr[i].X_fdb[0]-0.7f, lqr[i].X_fdb[0]+0.7f);
		//支持力解算
		float L0_array[3] = {vmc[i].L_fdb, vmc[i].V_fdb.e.vy0_fdb, vmc[i].Acc_fdb.L0_ddot};
		float theta_array[3] = {lqr[i].X_fdb[2], lqr[i].X_fdb[3], lqr[i].dot_x4};
		wlr.side[i].Fn_fdb = wlr_fn_calc(wlr.az_fdb, vmc[i].F_fdb.e.Fy_fdb, vmc[i].F_fdb.e.T0_fdb, L0_array, theta_array);
        kal_fn[i].measured_vector[0] = wlr.side[i].Fn_fdb;
        kalman_filter_update(&kal_fn[i]);
        wlr.side[i].Fn_kal = kal_fn[i].filter_vector[0];
		//离地检测
        if (wlr.high_flag != 2) {
            if(wlr.side[i].Fn_kal < 20.0f)
                wlr.side[i].fly_cnt++;
            else if(wlr.side[i].fly_cnt > 0)
                wlr.side[i].fly_cnt-=2;
            if(wlr.side[i].fly_cnt > 30) {
                wlr.side[i].fly_cnt = 30;
                wlr.side[i].fly_flag = 1;
            } else if(wlr.side[i].fly_cnt == 0)
                wlr.side[i].fly_flag = 0;
        } else {
            wlr.side[i].fly_flag = 0;
            wlr.side[i].fly_cnt = 0;
        }
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
	} else if (wlr.high_flag == 2) {
        wlr.high_set = LegLengthHigh2;
        pid_leg_length[0].kp = 800;
        pid_leg_length[1].kp = 800;
        pid_leg_length[0].kd = 30000;
        pid_leg_length[1].kd = 30000;
    } else if (wlr.high_flag == 1) { //长腿长
        pid_leg_length[0].kp = 500;
        pid_leg_length[1].kp = 500;
        pid_leg_length[0].kd = 10000;
        pid_leg_length[1].kd = 10000;
        if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {//腾空
            wlr.high_set = LegLengthHightFly;
        } else {
            wlr.high_set = LegLengthHigh;
        }
    } else { //正常腿长
        pid_leg_length[0].kp = 500;
        pid_leg_length[1].kp = 500;
        pid_leg_length[0].kd = 10000;
        pid_leg_length[1].kd = 10000;
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
//            k_array_fit1(K_Array_test, vmc[0].L_fdb);
//            aMartix_Cover(lqr[0].K, (float*)K_Array_test, 2, 6);
//            k_array_fit1(K_Array_test, vmc[1].L_fdb);
//            aMartix_Cover(lqr[1].K, (float*)K_Array_test, 2, 6);
            k_array_fit2(K_Array_test, vmc[0].L_fdb, vmc[0].q_fdb[0]/PI*180);
            aMartix_Cover(lqr[0].K, (float*)K_Array_test, 2, 6);
            k_array_fit2(K_Array_test, vmc[1].L_fdb, vmc[1].q_fdb[0]/PI*180);
            aMartix_Cover(lqr[1].K, (float*)K_Array_test, 2, 6);
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
        data_limit(&wlr.side[i].Tw, -4.5f, 4.5f);
	}
}

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
const float BodyWidth = 0.42f;//两轮间距
const float WheelRadius = 0.06f;//轮子半径
const float LegLengthMax = 0.30f, LegLengthMin = 0.10f;

const float LegLengthJump1 = 0.15f;//压腿
const float LegLengthJump2 = 0.35f;//蹬腿
const float LegLengthJump3 = 0.24f;//收腿
const float LegLengthJump4 = 0.22f;//落地

const float LegLengthHightFly = 0.28f;//长腿腿长腾空 0.28
const float LegLengthFly = 0.23f;//正常腿长腾空
const float LegLengthHigh = 0.23f;//长腿 0.23
const float LegLengthNormal = 0.18f;//正常

float x3_balance_zero = 0.08f, x5_balance_zero = 0.02f;//腿摆角角度偏置 机体俯仰角度偏置
float x3_fight_zero = 0.06f;
//								位移  速度	角度	角速度  角度	角速度
float K_Array_Wheel[2][6] =		{{60, 30, 80, 8, 300, 10}, 
                                { -0, -0.7, -8, -1, 3, 2}};
float K_Array_Leg[2][6] =		{{6.08535, 12.6636, 49.1418, 7.57203, 54.8073, 12.7387}, {-0.793527, -1.75976, -13.1544, -1.78789, 4.8796, 1.66868}};
float K_Array_Fly[2][6] =		{{0, 0, 80, 10, 0, 0}, { 0, 0, 0, 0, 0, 0}};
float K_Array_Prone[2][6] =     {{0, 0, 0, 0, 0, 0}, {0, 20, 0, 0, 0, 2}};

float K_Array_test[2][6];
float K_Array_List[12][4] = 
//{{6.33331,37.1637,-173.267,195.695},{13.375,62.6818,-306.463,349.224},{14.3574,523.661,-1613.6,1529.67},{4.56946,50.9723,-144.236,129.666},{32.2546,24.3717,604.921,-950.859},{9.25451,-6.78675,173.053,-249.371},{-0.907961,1.91211,-9.19263,11.1016},{-2.00572,5.03316,-22.1698,25.8241},{-4.53142,-42.3224,36.5561,-9.40343},{-0.980517,-1.26602,-10.414,10.7121},{13.699,-30.8227,43.8683,-31.5937},{2.94063,-5.26768,6.22139,-3.89125}};    
//{{7.7883,46.688,-215.454,242.688},{15.0481,71.2191,-344.839,392.009},{16.0739,537.105,-1644.75,1553.31},{4.70938,54.5476,-149.965,133.102},{29.3234,27.6173,614.168,-967.687},{8.55949,-6.85075,177.838,-255.908},{-1.11401,2.46758,-11.7282,14.1079},{-2.26846,5.99326,-26.0743,30.2096},{-4.74503,-41.8374,32.0378,-4.75168},{-1.0219,-1.15369,-11.5588,11.717},{14.0627,-32.6827,49.9302,-38.2469},{3.02982,-5.67978,7.63334,-5.45232}};
{{7.00681,41.1423,-185.674,205.977},{15.4159,74.7049,-349.699,390.408},{14.0718,685.814,-2051.,1902.38},{5.08813,65.8114,-175.561,152.294},{18.566,-19.1885,902.991,-1325.06},{6.46489,-20.1441,248.876,-340.781},{-0.730768,2.12264,-9.87338,11.7523},{-1.68527,5.33373,-23.4529,27.2529},{-4.39641,-36.0867,18.7682,8.33676},{-0.929593,-0.55816,-11.7229,12.0788},{12.6942,-29.997,50.1496,-42.3926},{2.70374,-5.18173,7.91857,-6.59155}};

float K_Array_List_lq[12][4][4] = 
{{{5.22798,0.120698,-0.000840624,1.30043e-6},{-27.8965,0.476016,-0.00280186,0},{-63.4158,0.0760932,0,0},{116.423,0,0,0}},
{{8.71734,0.29992,-0.00203995,2.86743e-6},{-71.0269,0.76054,-0.00457909,0},{-15.2343,0.168384,0,0},{71.7744,0,0,0}},
{{-22.99,1.69686,-0.0119447,0.0000182525},{-99.7898,4.6908,-0.0268931,0},{-761.51,0.472008,0,0},{1096.43,0,0,0}},
{{5.83449,0.0891018,-0.000731637,1.77972e-6},{-48.1618,0.816036,-0.00471951,0},{31.5375,0.0972437,0,0},{-53.7863,0,0,0}},
{{75.0024,-1.07961,0.00679659,-6.20881e-6},{144.783,1.05686,-0.00505925,0},{-328.529,-0.3575,0,0},{217.495,0,0,0}},
{{15.5352,-0.311456,0.00198484,-2.0315e-6},{75.1284,0.202374,-0.000786395,0},{-187.928,-0.148803,0,0},{181.46,0,0,0}},
{{-1.33258,0.0159886,-0.000100859,9.22646e-8},{-1.41807,-0.0127862,0.0000605034,0},{1.05816,0.00470869,0,0},{2.29189,0,0,0}},
{{-2.72025,0.0238647,-0.000148877,1.18866e-7},{-0.375646,0.0154076,-0.0000903077,0},{-14.2395,0.00189886,0,0},{26.0355,0,0,0}},
{{-4.91795,-0.0760693,0.00086145,-3.06738e-6},{54.8722,-2.12031,0.0117915,0},{77.9843,-0.0406164,0,0},{-84.3636,0,0,0}},
{{-2.00836,0.00600775,-0.0000267199,-5.19685e-8},{-0.986647,0.118254,-0.000643844,0},{-33.8391,-0.00886079,0,0},{46.0011,0,0,0}},
{{4.86264,0.0699955,-0.00050198,8.93186e-7},{-33.212,0.353139,-0.00210746,0},{12.5199,0.0682708,0,0},{0.422155,0,0,0}},
{{2.93993,0.00177889,-0.0000390133,2.43053e-7},{-14.3874,0.104903,-0.000639501,0},{20.4728,0.026013,0,0},{-23.4521,0,0,0}}};

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
	wlr.q0_set = PI / 2;
	
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
//        kal_v[i].Q_data[0] = 0.0003f; kal_v[i].Q_data[3] = 0.01f;
//        kal_v[i].R_data[0] = 1.0f; kal_v[i].R_data[3] = 10.0f;
        kal_v[i].Q_data[0] = 0.01f; kal_v[i].Q_data[3] = 0.01f;
        kal_v[i].R_data[0] = 1.0f; kal_v[i].R_data[3] = 10.0f;
		//PID参数初始化
		pid_init(&pid_leg_length[i], NONE, 500, 0.0f, 20000, 10, 20);//i 2.5f
		pid_init(&pid_leg_length_fast[i], NONE, 1000, 0, 10000, 30, 50);
	}
	//卡尔曼滤波器初始化

	//PID参数初始化
    pid_init(&pid_yaw, NONE, -7, 0, 0, 0, 10);
	pid_init(&pid_wz, NONE, 2.0f, 0, 7.0f, 0, 2.5f);
	pid_init(&pid_q0, NONE, 60, 0, 100, 0, 10);//与LQR的虚拟腿摆角控制拮抗 60 0 120
	pid_init(&pid_roll, NONE, 500, 0, 3000, 0, 30);//与VMC的腿长控制协同  1000 0 3500
//    pid_init(&pid_roll, NONE, 0, 0, 0, 0, 0);
//    pid_init(&pid_wx, NONE, 13, 0, 0, 0, 10);
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
    float v_fdb = (kal_v[0].filter_vector[0]+kal_v[1].filter_vector[0])/2.0f;
    if (fabs(v_fdb) > fabs(wlr.v_set))//加强超速控制
        wlr.v_set = data_fusion(wlr.v_set, 0, fabs(v_fdb - wlr.v_set));
    if (chassis.mode == CHASSIS_MODE_REMOTER_ROTATE ||
        chassis.mode == CHASSIS_MODE_KEYBOARD_ROTATE) {
        wlr.wz_set = 9;                                                 //小陀螺时固定转速 减小前进
        wlr.v_set = -0.5f * wlr.v_set;
    } else {
        wlr.wz_set = pid_calc(&pid_yaw, wlr.yaw_fdb + wlr.yaw_err, wlr.yaw_fdb);
        wlr.v_set = data_fusion(wlr.v_set, 0, fabs(wlr.wz_fdb)/6.0f);   //限制前进速度 旋转过快
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
        //定点控制第一方案
//        lqr[i].X_fdb[0] += (lqr[i].X_fdb[1] + lqr[i].last_x2) / 2 * CHASSIS_PERIOD_DU * 0.04f;//使用梯形积分速度求位移
        //定点控制第二方案
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
            //定点控制第一方案	
//            lqr[i].X_fdb[0] = 0;
            //定点控制第二方案
            lqr[i].X_ref[0] = lqr[i].X_fdb[0];
        }
		//支持力解算
		float L0_array[3] = {vmc[i].L_fdb, vmc[i].V_fdb.e.vy0_fdb, vmc[i].Acc_fdb.L0_ddot};
		float theta_array[3] = {lqr[i].X_fdb[2], lqr[i].X_fdb[3], lqr[i].dot_x4};
		wlr.side[i].Fn_fdb = wlr_fn_calc(wlr.az_fdb, vmc[i].F_fdb.e.Fy_fdb, vmc[i].F_fdb.e.T0_fdb, L0_array, theta_array);
        kal_fn[i].measured_vector[0] = wlr.side[i].Fn_fdb;
        kalman_filter_update(&kal_fn[i]);
        wlr.side[i].Fn_kal = kal_fn[i].filter_vector[0];
		//离地检测
        if(wlr.side[i].Fn_kal < 30.0f)
            wlr.side[i].fly_cnt++;
        else if(wlr.side[i].fly_cnt != 0)
            wlr.side[i].fly_cnt--;
        if(wlr.side[i].fly_cnt > 30) {
            wlr.side[i].fly_cnt = 30;
            wlr.side[i].fly_flag = 1;
        } else if(wlr.side[i].fly_cnt == 0)
        wlr.side[i].fly_flag = 0;
	}
	//高度选择 跳跃状态改变
	if (wlr.jump_flag == 1) {//跳跃起跳状态 先压腿
		wlr.high_set = LegLengthJump1;
		if(vmc[0].L_fdb < LegLengthJump1 && vmc[1].L_fdb < LegLengthJump1)
			wlr.jump_flag = 2;
	} else if (wlr.jump_flag == 2) {//起跳 弹腿
		wlr.high_set = LegLengthJump2;
		if(vmc[0].L_fdb > LegLengthJump2 && vmc[1].L_fdb > LegLengthJump2)
			wlr.jump_flag = 3;
	} else if (wlr.jump_flag == 3) {//收腿
		wlr.high_set = LegLengthJump3;
		if (vmc[0].L_fdb < LegLengthJump3 && vmc[1].L_fdb < LegLengthJump3 && !wlr.side[0].fly_flag && !wlr.side[1].fly_flag)
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
            wlr.high_set = LegLengthFly;
        } else {
            wlr.high_set = LegLengthHigh;
        }
    } else { //正常腿长
        wlr.high_set = LegLengthNormal;
    }
    //限制 旋转压腿长用于解决侧翻 腿摆压腿长用于解决膝关节可能卡地面
    if (wlr.high_flag) {
        wlr.high_set = data_fusion(wlr.high_set, 0.90f * wlr.high_set, fabs(wlr.wz_set/6.0f));
    }
    wlr.high_set = data_fusion(wlr.high_set, LegLengthMin, (fabs(lqr[0].X_fdb[2])+fabs(lqr[1].X_fdb[2]))/2.0f/2.0f);
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
        wlr.q0_offs   = pid_calc(&pid_q0, vmc[0].q_fdb[0], vmc[1].q_fdb[0]);//双腿摆角同步控制
    }
	wlr.roll_offs = pid_calc(&pid_roll, wlr.roll_set, wlr.roll_fdb);
//    wlr.wx_set = pid_calc(&pid_roll, wlr.roll_set, wlr.roll_fdb);
//    wlr.roll_offs = pid_calc(&pid_wx, wlr.wx_set, wlr.wx_fdb);
    if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {//腾空
        wlr.wz_offs   = 0;
    } else {
        wlr.wz_offs   = pid_calc(&pid_wz, wlr.wz_fdb, wlr.wz_set);//Yaw控制
    }
	//两侧轮腿分别独立控制
	for (int i = 0; i < 2; i++) {
		//LQR输入控制值 计算得出轮子与腿的力矩
        //速度融合 在快要撞限位时减小速度控制
        if (wlr.side[i].q1 > 3.6f)
            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (wlr.side[i].q1 - 3.6f)/0.3f);
        else if (wlr.side[i].q4 < -0.3f)
            lqr[i].X_ref[1] = data_fusion(twm.v_ref[i], lqr[i].X_fdb[1], (-0.3f - wlr.side[i].q4)/0.3f);
        else
            lqr[i].X_ref[1] = twm.v_ref[i];
        //摔倒保护 在快要摔倒时减小速度输入
		aMartix_Add(1, lqr[i].X_ref, -1, lqr[i].X_fdb, lqr[i].X_diff, 6, 1);
		aMartix_Mul(lqr[i].K, lqr[i].X_diff, lqr[i].U_ref, 2, 6, 1);
		//腿部虚拟力控制
		if (wlr.jump_flag == 2 || wlr.jump_flag == 4)						//跳跃蹬腿阶段 响应要大
			wlr.side[i].Fy = pid_calc(&pid_leg_length_fast[i], tlm.l_ref[i], vmc[i].L_fdb)\
								+ mb * GRAVITY / 2 + WLR_SIGN(i) * wlr.roll_offs;
		else if (wlr.jump_flag == 3)                                        //跳跃收腿阶段 响应要大
			wlr.side[i].Fy = pid_calc(&pid_leg_length_fast[i], tlm.l_ref[i], vmc[i].L_fdb);
		else if (wlr.side[0].fly_flag && wlr.side[1].fly_flag) {            //浮空收腿 响应不用那么大
			wlr.side[i].Fy = pid_calc(&pid_leg_length_fast[i], tlm.l_ref[i], vmc[i].L_fdb);
		} else																//常态 跳跃压腿阶段 跳跃落地阶段
			wlr.side[i].Fy = pid_calc(&pid_leg_length[i], tlm.l_ref[i], vmc[i].L_fdb)\
                                 + 27 + WLR_SIGN(i) * wlr.roll_offs;
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

#include "gimbal_task.h"
#include "mode_switch_task.h"
#include "control_def.h"
#include "drv_dji_motor.h"
#include "prot_vision.h"
#include "prot_dr16.h"
#include "prot_imu.h"
#include "math_lib.h"
#include "cmsis_os.h"
#include "string.h"

gimbal_scale_t gimbal_scale = {
    .ecd_remote = 1,
    .ecd_keyboard = 1,
    .angle_remote = 0.00002f,
    .angle_keyboard = 0.00006f
};
gimbal_t gimbal;
uint8_t test_gimbal_vision_mode = 0;//遥控0 视觉1

static void gimbal_init(void)
{
    memset(&gimbal, 0, sizeof(gimbal_t));
    pid_init(&gimbal.yaw_angle.pid, 30, 0, 400, 0, 10);
    pid_init(&gimbal.yaw_spd.pid, 0.5f, 0.006f, 0, 0.6f, 1.3f);
//    pid_init(&gimbal.yaw_ecd.pid, 0, 0, 0, 0, 0);

    pid_init(&gimbal.pit_angle.pid, 30, 0, 500, 0, 15);
    pid_init(&gimbal.pit_spd.pid, -0.2f, -0.003f, 0, 0.6f, 1.3f);
//    pid_init(&gimbal.pit_ecd.pid, 0, 0, 0, 0, 0);
}

static void gimbal_pid_calc(void)
{
    //位置环反馈 陀螺仪
    //速度环反馈 陀螺仪
    data_limit(&gimbal.pit_angle.ref, -0.7f, 0.3f);
    gimbal.pit_angle.fdb = gimbal_imu.pit;
    gimbal.pit_spd.ref = pid_calc(&gimbal.pit_angle.pid, gimbal.pit_angle.ref, gimbal.pit_angle.fdb);
    gimbal.pit_spd.fdb = gimbal_imu.wy;
    gimbal.pit_output = pid_calc(&gimbal.pit_spd.pid, gimbal.pit_spd.ref, gimbal.pit_spd.fdb);

    if (gimbal.yaw_angle.ref < 0) {
        gimbal.yaw_angle.ref += 2 * PI;
    } else if (gimbal.yaw_angle.ref > 2 * PI) {
        gimbal.yaw_angle.ref -= 2 * PI;
    }
    gimbal.yaw_angle.fdb = gimbal_imu.yaw;
    float yaw_err = circle_error(gimbal.yaw_angle.ref, gimbal.yaw_angle.fdb, 2*PI);
    gimbal.yaw_spd.ref = pid_calc(&gimbal.yaw_angle.pid, gimbal.yaw_angle.fdb + yaw_err, gimbal.yaw_angle.fdb);
    gimbal.yaw_spd.fdb = gimbal_imu.wz;
    gimbal.yaw_output = -0.07f + pid_calc(&gimbal.yaw_spd.pid, gimbal.yaw_spd.ref, gimbal.yaw_spd.fdb);
    //位置环反馈 电机编码器
    //速度环反馈 陀螺仪
}

static void gimbal_data_output(void)
{
    dji_motor_set_torque(&pit_motor, gimbal.pit_output);
    dji_motor_set_torque(&yaw_motor, gimbal.yaw_output);
}

static void gimbal_get_vision_data(void)
{
    switch (vision.aim_status) {
        case AIMING: {//识别到目标
            if (vision.new_frame_flag) {
                vision.new_frame_flag = 0;
                gimbal.pit_angle.ref = vision.target_pit_angle;
                gimbal.yaw_angle.ref = vision.target_yaw_angle;
            }
            break;
        }
        case FIRST_LOST: {//首次丢失
            vision.aim_status = UNAIMING;
            gimbal.pit_angle.ref = gimbal.pit_angle.fdb;
            gimbal.yaw_angle.ref = gimbal.yaw_angle.fdb;
            break;
        }
        case UNAIMING: {//未识别到目标
            if (rc.sw1 == RC_MI) {
//                gimbal.pit_ecd.ref -= rc.ch2 * gimbal_scale.ecd_remote;
                gimbal.pit_angle.ref -= rc.ch2 * gimbal_scale.angle_remote;
//                gimbal.yaw_ecd.ref -= rc.ch1 * gimbal_scale.ecd_remote;
                gimbal.yaw_angle.ref -= rc.ch1 * gimbal_scale.angle_remote;
            } else if (rc.sw1 == RC_DN) {
                gimbal.pit_angle.ref += rc.mouse.y * gimbal_scale.angle_keyboard;
                gimbal.yaw_angle.ref -= rc.mouse.x * gimbal_scale.angle_keyboard;
            }
            break;
        }
        default: break;
    }
}

void gimbal_task(void const *argu)
{
    uint32_t thread_wake_time = osKernelSysTick();
    gimbal_init();
    for(;;) {
//        taskENTER_CRITICAL();
        switch (ctrl_mode) {
            case PROTECT_MODE: {
                gimbal.yaw_angle.ref = gimbal_imu.yaw;
                gimbal.pit_angle.ref = 0;
//                gimbal.pit_ecd.ref = GIMBAL_PIT_CENTER_OFFSET;
//                gimbal.yaw_ecd.ref = GIMBAL_YAW_CENTER_OFFSET;
                gimbal.pit_output = 0;
                gimbal.yaw_output = 0;
                break;
            }
            case REMOTER_MODE: {
//                gimbal.pit_ecd.ref -= rc.ch2 * gimbal_scale.ecd_remote;
                gimbal.pit_angle.ref -= rc.ch2 * gimbal_scale.angle_remote;
//                gimbal.yaw_ecd.ref -= rc.ch1 * gimbal_scale.ecd_remote;
                gimbal.yaw_angle.ref -= rc.ch1 * gimbal_scale.angle_remote;
                gimbal_pid_calc();
                break;
            }
            case KEYBOARD_MODE: {
                if (kb_status[KEY_SHOOT_HOUSE]) { //补弹
                    gimbal.pit_angle.ref = 0;
                    gimbal.yaw_angle.ref -= rc.mouse.x * gimbal_scale.angle_keyboard * 0.5f;
                } else {
                    //一键调头
                    if(key_scan_clear(KEY_GIMBAL_TURN_R)) {
                        gimbal.yaw_angle.ref -= PI/2;
                    } 
                    else if (key_scan_clear(KEY_GIMBAL_TURN_L)) {
                        gimbal.yaw_angle.ref += PI/2;
                    }
                    gimbal.pit_angle.ref += rc.mouse.y * gimbal_scale.angle_keyboard * 0.5f;
                    gimbal.yaw_angle.ref -= rc.mouse.x * gimbal_scale.angle_keyboard;
                }
                gimbal_pid_calc();
                break;
            }
            case VISION_MODE: { //视觉模式 
                gimbal_get_vision_data();
                gimbal_pid_calc();
                break;
            }
            default:break;
        }
        gimbal_data_output();
//        taskEXIT_CRITICAL();
        osDelayUntil(&thread_wake_time, 2);
    }
}

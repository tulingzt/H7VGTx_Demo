#include "prot_judge.h"
#include "cmsis_os.h"
#include "string.h"
#include "crc.h"

//机器人阵营
typedef enum
{
    BLUE,
    RED
} color_e;

UART_HandleTypeDef *judge_huart;
frame_header_t frame_header;

//机器人接收的数据
game_status_t                       game_status;                    //比赛状态数据，固定以1Hz频率发送
game_result_t                       game_result;                    //比赛结果数据，比赛结束触发发送
game_robot_HP_t                     game_robot_HP;                  //机器人血量数据，固定以3Hz频率发送
event_data_t                        event_data;                     //场地事件数据，固定以1Hz频率发送
ext_supply_projectile_action_t      ext_supply_projectile_action;   //补给站动作标识数据，补给站弹丸释放时触发发送
referee_warning_t                   referee_warning;                //裁判警告数据，己方判罚/判负时触发发送，其余时间以1Hz频率发送
dart_info_t                         dart_info;                      //飞镖发射相关数据，固定以1Hz频率发送
robot_status_t                      robot_status;                   //机器人性能体系数据，固定以10Hz频率发送
power_heat_data_t                   power_heat_data;                //实时底盘功率和枪口热量数据，固定以50Hz频率发送
robot_pos_t                         robot_pos;                      //机器人位置数据，固定以1Hz频率发送
buff_t                              buff;                           //机器人增益数据，固定以3Hz频率发送
hurt_data_t                         hurt_data;                      //伤害状态数据，伤害发生后发送
shoot_data_t                        shoot_data;                     //实时射击数据，弹丸发射后发送
projectile_allowance_t              projectile_allowance;           //允许发弹量，固定以10Hz频率发送
rfid_status_t                       rfid_status;                    //机器人RFID模块状态，固定以3Hz频率发送
////空中独有的接收
//air_support_data_t                  air_support_data;               //空中支援时间数据，固定以1Hz频率发送
////飞镖独有的接收
//dart_client_cmd_t                   dart_client_cmd;                //飞镖选手端指令数据，固定以3Hz频率发送
////哨兵独有的接收
//ground_robot_position_t             ground_robot_position;          //地面机器人位置数据，固定以1Hz频率发送
//sentry_info_t                       sentry_info;                    //哨兵自主决策信息同步，固定以1Hz频率发送
////雷达独有的接收
//radar_mark_data_t                   radar_mark_data;                //雷达标记进度数据，固定以1Hz频率发送
//radar_info_t                        radar_info;                     //雷达自主决策信息同步，固定以1Hz频率发送

//机器人交互/ui
robot_interaction_data_t            robot_interaction_data;         //机器人交互数据，发送方触发发送，频率上限为10Hz
//小地图交互接收
map_command_t                       map_command;                    //选手端小地图交互数据，选手端触发发送
map_robot_data_t                    map_robot_data;                 //选手端小地图接收雷达数据，频率上限为10Hz
//小地图交互发送
//map_data_t                          map_data;                       //选手端小地图接收哨兵数据，频率上限为1Hz
//custom_info_t                       custom_info;                    //选手端小地图接收机器人数据，频率上限为3Hz
////自定义控制器接收
//custom_robot_data_t                 custom_robot_data;              //自定义控制器与机器人交互数据，发送方触发发送，频率上限为30Hz        图传链路
////自定义控制器发送
//custom_client_data_t                custom_client_data;             //自定义控制器与选手端交互数据，发送方触发发送，频率上限为30Hz
//图传链路遥控器
remote_control_t                    remote_control;                 //键鼠遥控数据，固定30Hz频率发送                                      图传链路

/*
 * @brief     读取裁判数据函数，串口中断函数中直接调用进行读取
 * @param[in] data: 数据指针
 * @retval    数据正常返回0，异常返回1
 */
uint8_t judge_get_data(uint8_t *data)
{
    uint8_t result = 1;
    uint16_t data_length;
    int cmd_id;
    //写入帧头
//    memcpy(&frame_header, data, 5);
    frame_header.SOF = data[0];
    frame_header.data_length = (data[2] << 8 | data[1]);
    frame_header.seq = data[3];
    frame_header.CRC8 = data[4];
    //判断帧头数据是否为0xA5
    if (frame_header.SOF == 0xA5) {
        //帧头CRC8校验
        if (crc8_verify_checksum(data, 5)) {
            //统计一帧数据长度,用于CRC16校验
            data_length = frame_header.data_length + 5 + 2 + 2;
            //帧尾CRC16校验
            if (crc16_verify_checksum(data, data_length)) {
                result = 0;
                cmd_id = (data[6] << 8 | data[5]);
                switch (cmd_id) {
                        case ID_game_status                  : memcpy(&game_status                 , (data + 7), LEN_game_status                 );break;
                        case ID_game_result                  : memcpy(&game_result                 , (data + 7), LEN_game_result                 );break;
                        case ID_game_robot_HP                : memcpy(&game_robot_HP               , (data + 7), LEN_game_robot_HP               );break;
                        case ID_event_data                   : memcpy(&event_data                  , (data + 7), LEN_event_data                  );break;
                        case ID_ext_supply_projectile_action : memcpy(&ext_supply_projectile_action, (data + 7), LEN_ext_supply_projectile_action);break;
                        case ID_referee_warning              : memcpy(&referee_warning             , (data + 7), LEN_referee_warning             );break;
                        case ID_dart_info                    : memcpy(&dart_info                   , (data + 7), LEN_dart_info                   );break;
                        case ID_robot_status                 : memcpy(&robot_status                , (data + 7), LEN_robot_status                );break;
                        case ID_power_heat_data              : memcpy(&power_heat_data             , (data + 7), LEN_power_heat_data             );break;
                        case ID_robot_pos                    : memcpy(&robot_pos                   , (data + 7), LEN_robot_pos                   );break;
                        case ID_buff                         : memcpy(&buff                        , (data + 7), LEN_buff                        );break;
//                        case ID_air_support_data             : memcpy(&air_support_data            , (data + 7), LEN_air_support_data            );break;
                        case ID_hurt_data                    : memcpy(&hurt_data                   , (data + 7), LEN_hurt_data                   );break;
                        case ID_shoot_data                   : memcpy(&shoot_data                  , (data + 7), LEN_shoot_data                  );break;
                        case ID_projectile_allowance         : memcpy(&projectile_allowance        , (data + 7), LEN_projectile_allowance        );break;
                        case ID_rfid_status                  : memcpy(&rfid_status                 , (data + 7), LEN_rfid_status                 );break;
//                        case ID_dart_client_cmd              : memcpy(&dart_client_cmd             , (data + 7), LEN_dart_client_cmd             );break;
//                        case ID_ground_robot_position        : memcpy(&ground_robot_position       , (data + 7), LEN_ground_robot_position       );break;
//                        case ID_radar_mark_data              : memcpy(&radar_mark_data             , (data + 7), LEN_radar_mark_data             );break;
//                        case ID_sentry_info                  : memcpy(&sentry_info                 , (data + 7), LEN_sentry_info                 );break;
//                        case ID_radar_info                   : memcpy(&radar_info                  , (data + 7), LEN_radar_info                  );break;
                        case ID_robot_interaction_data       : memcpy(&robot_interaction_data      , (data + 7), LEN_robot_interaction_data      );break;
//                        case ID_custom_robot_data            : memcpy(&custom_robot_data           , (data + 7), LEN_custom_robot_data           );break;
                        case ID_map_command                  : memcpy(&map_command                 , (data + 7), LEN_map_command                 );break;
                        case ID_remote_control               : memcpy(&remote_control              , (data + 7), LEN_remote_control              );break;
                        case ID_map_robot_data               : memcpy(&map_robot_data              , (data + 7), LEN_map_robot_data              );break;
//                        case ID_custom_client_data           : memcpy(&custom_client_data          , (data + 7), LEN_custom_client_data          );break;
//                        case ID_map_data                     : memcpy(&map_data                    , (data + 7), LEN_map_data                    );break;
//                        case ID_custom_info                  : memcpy(&custom_info                 , (data + 7), LEN_custom_info                 );break;
                }
            }
        }
        //首地址加帧长度,指向CRC16下一字节,用来判断是否为0xA5,用来判断一个数据包是否有多帧数据 
        if(*(data + sizeof(frame_header) + 2 + frame_header.data_length + 2) == 0xA5) {
            //如果一个数据包出现了多帧数据,则再次读取
            judge_get_data(data + sizeof(frame_header) + 2 + frame_header.data_length + 2);
        }
    }
    return result;
}

void judge_init(UART_HandleTypeDef *huart)
{
    judge_huart = huart;
}

/*
 * @brief     发送裁判数据
 * @param[in] frame_header: 帧头结构体
 * @param[in] cmd_id      : 帧id
 * @param[in] data        : 数据指针
 * @retval    void
 */
uint8_t send_data[150];
static void judge_send_data(frame_header_t *frame_header, cmd_id_e cmd_id, uint8_t *data)
{
    //如果串口对应的DMA还未发送完就等待，防止变量UartTxBuf被CPU和DMA同时使用。
    //发送裁判系统数据不用急
    while (judge_huart->gState != HAL_UART_STATE_READY) {
        osDelay(1);
    }
    memcpy(send_data, frame_header, 5);//写入帧头数据
    send_data[5] = cmd_id >> 8;//写入id
    send_data[6] = cmd_id & 0x0f;
    memcpy(send_data + 7, data, frame_header->data_length);
    crc16_set_checksum(send_data, frame_header->data_length + 9);
    HAL_UART_Transmit_DMA(judge_huart, send_data, frame_header->data_length + 9);
}

/*
 * @brief     发送机器人交互数据
 * @param[in] rcv_id: 接收方id
 * @param[in] id    : 子数据id
 * @param[in] len   : 子数据长度
 * @param[in] data  : 子数据指针
 * @retval    void
 */
static void robot_interaction_send_data(uint16_t rcv_id, robot_interaction_id_e id, robot_interaction_len_e len, uint8_t *data)
{
    frame_header_t frame_header;
    robot_interaction_data_t robot_interaction_data;
    frame_header.SOF = 0xA5;
    frame_header.data_length = len + 6;
    frame_header.seq = 0;
    crc8_set_checksum((uint8_t *)&frame_header, 5);
    robot_interaction_data.data_cmd_id = id;
    //发送接收id设置
    robot_interaction_data.sender_id = robot_status.robot_id;
    robot_interaction_data.receiver_id = rcv_id;
    
    memcpy((uint8_t *)&robot_interaction_data, data, len);
    
    judge_send_data(&frame_header, 0x0301, (uint8_t *)&robot_interaction_data);
}



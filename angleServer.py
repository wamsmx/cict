import json
import requests
import websocket
from requests_ntlm import HttpNtlmAuth as wauth
import time
import math
#define WebSocket-Client
ws = websocket.WebSocket(sslopt={"check_hostname": False})
#connect
ws.connect("ws://localhost:8000/ws/polData/")

#Routine to get the angle and calculate the relatie angle value
def correct_angle(value):
    value = round(value, 6)
    if value > 360:
        value = value - 360
    elif value < 0:
        value = value + 360
    else:
        value = value
    value = value*math.pi/180
    return round(value, 6)
names_dict = {'260': 'Tampico', '250':'Ciudad Juárez', '175': 'Morelia', '11': 'Tuxtla', '63': 'Guadalajara', '205': 'Cdmx', '85': 'Monterrey', '58': 'Mazatlan','53': 'Chetumal','68': 'El Salvador'}
_order = [175, 11, 63, 205, 85, 58, 53, 68, 260, 250]

session=requests.Session()
session.auth=wauth('potencia2021','W@ms_project202x')
url="http://localhost:6152/historian/timeseriesdata/read/current/175,11,63,205,85,58,53,68,260,250,290/json"

#This function calls the url and get the latest angle measurements
def get_angle():
    r=session.get(url)
    source=(r.content)
    source=json.loads(source)
    source = source['TimeSeriesDataPoints']
    angle = [d['Value'] for d in source]
    angle = list(map(lambda x : correct_angle(x), angle))
    return angle

#This function separates the angle array
def get_each_angle():
    angles = get_angle()
    morelia = angles[0] 
    tuxtla = angles[1] 
    guadalajara = angles[3] 
    cdmx = angles[3]   
    monterrey = angles[4] 
    mazatlan = angles[5]  
    chetumal = angles[6] 
    salvador = angles[7] 
    tampico = angles[8]
    juarez = angles[9]
    mnzlo = angles[10]
    return morelia, tuxtla, guadalajara, cdmx, monterrey, mazatlan, chetumal, salvador, tampico, juarez, mnzlo

angle0 = get_angle()
#Calculate ref angle
ref_cdmx_0 = angle0[3]
ref_chiapas_0 = angle0[1]

relative_angle_0 = []

for i in range(len(angle0)):
    if i == 7:
        _ref_angle_0 = ref_chiapas_0
    else:
        _ref_angle_0 = ref_cdmx_0

    _relative_agle_0 = angle0[i] - _ref_angle_0
    relative_angle_0.append(_relative_agle_0)

#wait a second
time.sleep(1)
n = 0
while True:
    if n == 1200:
        n = 0
        angle0 = get_angle()
        #Calculate ref angle
        ref_cdmx_0 = angle0[3]
        ref_chiapas_0 = angle0[1]

        relative_angle_0 = []

        for i in range(len(angle0)):
            if i == 7:
                _ref_angle_0 = ref_chiapas_0
            else:
                _ref_angle_0 = ref_cdmx_0

            _relative_agle_0 = angle0[i] - _ref_angle_0
            relative_angle_0.append(_relative_agle_0)

    #get angles
    angle1 = get_angle()

    '''Step 1: unwrap the angle measurement
        The angle measurements are limited within (0, 2PI), we need to 'unwrap' the angle to make the angle values continuous for further steps.
        Note: start from i = 1'''
    unwrapped_angles = []
    for i in range(len(angle0)):
        _angle0 = angle0[i]
        _angle1 = angle1[i]
        _unwrapped =  _angle1  - round((( _angle1  -  _angle0)/2/math.pi), 0)*2*math.pi
        unwrapped_angles.append(_unwrapped)

    angle0 = unwrapped_angles

    '''Step 2: Calculate relative phase angle
        Note: start from i = 0'''
    ref_cdmx = unwrapped_angles[3]
    ref_chiapas = unwrapped_angles[1]
    relative_angle = []
    for i in range(len(unwrapped_angles)):
        if i == 7:
            _ref_angle = ref_chiapas
        else:
            _ref_angle = ref_cdmx

        _relative_agle = unwrapped_angles[i] - _ref_angle
        relative_angle.append(_relative_agle)

    '''Step 3: Self-normalize the relative phase angle
    Shift the relative phase angle to make the angle series starts from 0.
    Note: start from i = 1'''  
    self_norm_angle = []  
    for i in range(len(relative_angle)):
        _self_norm_angle = round((relative_angle[i] - relative_angle_0[i])*180/math.pi,1)
        if abs(_self_norm_angle) > 90:
            _self_norm_angle = 0
        self_norm_angle.append(_self_norm_angle)

    '''dict_rel_angle ={'mor-cdmx': self_norm_angle[0], 'tuxt-cdmx': self_norm_angle[1], 'guad-cdmx': self_norm_angle[2],
        'cdmx-cdmx': self_norm_angle[3], 'mtrey-cdmx': self_norm_angle[4], 'mztln-cdmx': self_norm_angle[5],
        'chtmal-cdmx': self_norm_angle[6], 'slvd-tuxt': self_norm_angle[7], 'tmpc-cdmx': self_norm_angle[8],
        'cdjrz-cdmx': self_norm_angle[9]}'''
    
    self_norm_angle = [self_norm_angle[9], self_norm_angle[7],self_norm_angle[1], self_norm_angle[6], self_norm_angle[0],
                        self_norm_angle[2], self_norm_angle[5], self_norm_angle[4],self_norm_angle[8], self_norm_angle[10]]

    ws.send(json.dumps({'value': self_norm_angle}))
    #Wait 1 sec before iterate again
    time.sleep(1)
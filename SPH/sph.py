from vpython import *
import pyautogui

# Poly6 kernel
# def poly6Kernel(r, h):
#     if r < 0 or h < r: 
#         return 0
#     return (4 * (h * h - r * r) ** 3 / (pi * (h ** 8))) 
def poly6Kernel(r, h):
    if r < 0 or h < r: 
        return 0
    return (315 / (64 * pi * h ** 9)) * (h * h - r * r) ** 3

#Spiky kernel 
# def grid_Spiky(r, h, dis):
#     if r == 0 or r > h :
#         return 0 
#     return (30 / (pi * h ** 6) * (h - r) ** 2 * h * dis / r)
def grid_Spiky(r, h, dis):
    if r <= 0 or r > h:
        return 0 
    return (45 / (pi * h ** 6) * (h - r) ** 2) * dis
# viscosity kernel
# def viscosityKernel(r, h):
#     if r < 0 or r > h : 
#         return 0 
#     return (40 / (pi * h ** 5) * (h - r))
def viscosityKernel(r, h):
    if r < 0 or r > h:
        return 0 
    return (45 / (pi * h ** 5)) * (h - r)

# 경계면 처리 함수
def boundary(obj, y_min, y_max, x_min, x_max, repul):
    if obj.pos.y < y_min:
        obj.pos.y = y_min
        obj.vel.y *= -repul
    
    if obj.pos.y > y_max:
        obj.pos.y = y_max
        obj.vel.y *= -repul
        
    if obj.pos.x < x_min:
        obj.pos.x = x_min
        obj.vel.x *= -repul
        
    if obj.pos.x > x_max:
        obj.pos.x = x_max
        obj.vel.x *= -repul

# 화면 설정
scene.autoscale = True 
scene.center = vec(50,50,0)

# 상수 초기화
column = 15
row = 20

numOfSphere = column * row 

radius = 2
kernel = radius * 2
side_space = 0.9 * kernel
high_space = 0.9 * kernel

# 수조 만들기 & 파티클 만들기
watertank = box( pos = vec(55,50,0), size = vec(110, 150, 10), color = color.cyan, opacity = 0.2)

particles = [] 
for i in range(0, row) :
    for j in range(0, column) : 
        particles.append( sphere(pos = vec(2 + j * side_space, high_space * i - 22, 0)  , radius = radius, color = color.white ) )
    
    
# 물리 성질 초기화
for p in particles: 
    p.mass = 12
    p.vel = vec(0,0,0)
    p.density = 0
    p.pressure = 0
    p.viscosity = vec(0,0,0)
    p.pressureForce = vec(0,0,0)
    p.force = vec(0,0,0)
    
    
m_kernel_h = radius * 3 
m_limit_velocity = 80

mu = 100 # 마찰 계수
gravity = -9.8 # 중력 가속도
repulsive = 0.3 # 경계면 처리 상수
m_restDenstity = particles[0].mass * poly6Kernel(0, m_kernel_h) # 고유 밀도 
k = 5000 # 압력 밀도 관련 함수

#경계설정
boundary_xmin = watertank.pos.x - watertank.size.x / 2 + radius
boundary_xmax = watertank.pos.x + watertank.size.x / 2 - radius
boundary_ymin = watertank.pos.y - watertank.size.y / 2 + radius
boundary_ymax = watertank.pos.y + watertank.size.y / 2 - radius



# 시간 설정 
t = 0 
dt = 0.03
frame_count = 0
while t < 100000:
    rate(100) 
    if( frame_count % 5 == 0):
        pyautogui.screenshot("./movie/frame_%04d.png" % (frame_count / 5))
    #VPython 장면만 캡처
    #scene.capture("./movie/frame_%04d.png" % frame_count)

    # 밀도 업데이트
    for i in range(numOfSphere):
        rSum = 0
        for j in range(numOfSphere):
            rdistance = mag(particles[i].pos - particles[j].pos)
            if rdistance < 0 or rdistance > m_kernel_h:  # 너무 가깝거나 멀면 탈락
                continue
            
            # 밀도 * kernel
            rSum += particles[j].mass * poly6Kernel(rdistance, m_kernel_h)
        
        particles[i].density = rSum

    # 밀도로 본인 압력 구하기
    for i in range(numOfSphere):
        particles[i].pressure = k * (particles[i].density - m_restDenstity)
        
    # 압력힘
    for i in range(numOfSphere):
        psum = [0.0, 0.0, 0.0]
        for j in range(numOfSphere):
            pdistance = mag(particles[i].pos - particles[j].pos)
            if pdistance < 0 or pdistance > m_kernel_h : 
                continue
            
            psum[0] += particles[j].mass * (particles[i].pressure  + particles[j].pressure) * 0.5 / particles[j].density * grid_Spiky(pdistance, m_kernel_h, particles[i].pos.x - particles[j].pos.x)
            psum[1] += particles[j].mass * (particles[i].pressure  + particles[j].pressure) * 0.5 / particles[j].density * grid_Spiky(pdistance, m_kernel_h, particles[i].pos.y - particles[j].pos.y)
            psum[2] += particles[j].mass * (particles[i].pressure  + particles[j].pressure) * 0.5 / particles[j].density * grid_Spiky(pdistance, m_kernel_h, particles[i].pos.z - particles[j].pos.z)
        
        particles[i].pressureForce.x = psum[0]
        particles[i].pressureForce.y = psum[1]
        particles[i].pressureForce.z = 0
        
    # 점성힘
    for i in range(numOfSphere) :
        vsum = [0.0, 0.0, 0.0]
        for j in range (numOfSphere): 
            vdistance = mag( particles[i].pos - particles[j].pos ) 
            if vdistance < 0 or vdistance > m_kernel_h : 
                continue
            visKernel = viscosityKernel(vdistance, m_kernel_h) 
            vsum[0] += particles[j].mass * ( particles[j].vel.x - particles[i].vel.x )  / particles[j].density * visKernel 
            vsum[1] += particles[j].mass * ( particles[j].vel.y - particles[i].vel.y )  / particles[j].density * visKernel 
            vsum[2] += particles[j].mass * ( particles[j].vel.z - particles[i].vel.z )  / particles[j].density * visKernel   
        
        particles[i].viscosity.x = mu * vsum[0]
        particles[i].viscosity.y = mu * vsum[1]
        particles[i].viscosity.z = 0 #vsum[0]
        
        
    # 외력 
    for i in range(numOfSphere): 
        particles[i].force.x = particles[i].pressureForce.x + particles[i].viscosity.x
        particles[i].force.y = particles[i].pressureForce.y + particles[i].viscosity.y + particles[i].mass * gravity 
        
        particles[i].force.z = 0
        
        
    # 속도 업데이트
    for i in range(numOfSphere):
        particles[i].vel.x = particles[i].vel.x + particles[i].force.x / particles[i].mass * dt
        particles[i].vel.y = particles[i].vel.y + particles[i].force.y / particles[i].mass * dt
        particles[i].vel.z = 0.0
        
        if particles[i].vel.x > m_limit_velocity :
            particles[i].vel.x = m_limit_velocity
        if particles[i].vel.y > m_limit_velocity :
            particles[i].vel.y = m_limit_velocity
        if particles[i].vel.x < -m_limit_velocity :
            particles[i].vel.x = -m_limit_velocity
        if particles[i].vel.y < -m_limit_velocity :
            particles[i].vel.y = -m_limit_velocity
    
        # 경계면 처리
        boundary(particles[i], boundary_ymin, boundary_ymax, boundary_xmin, boundary_xmax, repulsive) 

        particles[i].pos = particles[i].pos + particles[i].vel * dt
        
    t += dt
    frame_count += 1
    
   
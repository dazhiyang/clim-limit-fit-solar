import os
import pandas as pd
from pvlib.iotools import get_cams
import time
from tqdm import tqdm

# Settings 
# 设置
EMAIL = "yangdazhi.nus@gmail.com"
BASE_DIR = "/Volumes/Macintosh Research/Data/BSRN_QC_isolation/Data"
METADATA_FILE = os.path.join(BASE_DIR, "BSRN information.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "McClear")
SERVER = "api.soda-solardata.com"

# Ensure output directory exists 
# 确保输出目录存在
os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_yymm(x):
    """Parse YYMM string into a 4-digit year. / 将 YYMM 字符串解析为 4 位年份。"""
    s = str(x).zfill(4)
    yr = int(s[:2])
    # Assume 19xx for yr > 80, else 20xx 
    # 假设 yr > 80 为 19xx，否则为 20xx
    full_yr = 1900 + yr if yr > 80 else 2000 + yr
    return full_yr

def main():
    # Load metadata for coordinates 
    # 加载元数据以获取坐标
    df_meta = pd.read_csv(METADATA_FILE)
    df_meta['stn_lower'] = df_meta['stn'].str.lower()
    
    # Path to BSRN data 
    # BSRN 数据路径
    bsrn_dir = os.path.join(BASE_DIR, "BSRN")
    
    # Get all potential station directories 
    # 获取所有潜在的站点目录
    if not os.path.exists(bsrn_dir):
        print(f"Error: BSRN directory not found at {bsrn_dir}")
        return
        
    stations = sorted([d for d in os.listdir(bsrn_dir) if os.path.isdir(os.path.join(bsrn_dir, d))])

    for index, stn_id in enumerate(stations):
        # Match with metadata for coordinates 
        # 与元数据匹配以获取坐标
        meta_row = df_meta[df_meta['stn_lower'] == stn_id]
        if meta_row.empty:
            print(f"Warning: No metadata found for station {stn_id}. Skipping.")
            continue
            
        lat = meta_row.iloc[0]['Latitude']
        lon = meta_row.iloc[0]['Longitude']
        
        # Determine years from file names 
        # 从文件名确定年份
        stn_bsrn_path = os.path.join(bsrn_dir, stn_id)
        files = [f for f in os.listdir(stn_bsrn_path) if f.endswith('.dat.gz')]
        
        years_set = set()
        for f in files:
            # Pattern: <stn><mm><yy>.dat.gz (e.g., ale0110.dat.gz)
            # Extract yy from the filename (7 chars for .dat.gz, then 2 for yy)
            # 文件名模式：<站点><月><年>.dat.gz (例如，ale0110.dat.gz)
            # 提取 yy 从文件名（7 个字符用于 .dat.gz，然后 2 个字符用于 yy）
            try:
                # Use f[-9:-7] to get the yy part
                # 使用 f[-9:-7] 获取 yy 部分
                yy_str = f[-9:-7]
                yy = int(yy_str)
                full_yr = 1900 + yy if yy > 80 else 2000 + yy
                years_set.add(full_yr)
            except (ValueError, IndexError):
                continue
        
        years = sorted(list(years_set))
        
        if not years:
            print(f"No valid years found for station {stn_id}. Skipping.")
            continue
            
        print(f"\nProcessing Station: {stn_id} ({index+1}/{len(stations)})")
        
        # Station folder structure (lowercase) 
        # 站点文件夹结构（小写）
        stn_dir = os.path.join(OUTPUT_DIR, stn_id)
        os.makedirs(stn_dir, exist_ok=True)
        
        for yr in tqdm(years, desc=f"Downloading {stn_id}"):
            out_file = os.path.join(stn_dir, f"{stn_id}_{yr}.csv")
            
            if os.path.exists(out_file):
                continue
            
            start_date = f"{yr}-01-01"
            end_date = f"{yr}-12-31"
            
            # Retry logic 
            # 重试逻辑
            max_retries = 5
            retry_delay = 5  
            
            for attempt in range(max_retries):
                try:
                    # get_cams returns a tuple (data, metadata) 
                    # time_step='1min' for 1-minute resolution 
                    # integrated=False for instantaneous values 
                    # get_cams 返回一个元组 (数据, 元数据)
                    # time_step='1min' 表示 1 分钟分辨率
                    # integrated=False 表示瞬时值
                    data, meta = get_cams(
                        latitude=lat,
                        longitude=lon,
                        start=start_date,
                        end=end_date,
                        email=EMAIL,
                        time_step='1min',
                        integrated=False,
                        url=SERVER
                    )
                    
                    if not data.empty:
                        if 'Observation period' in data.columns:
                            # Extract "END" time from "START/END" format (split by '/' and take second part)
                            # 从 "START/END" 格式中提取 "结束" 时间（按 '/' 分割并取第二部分）
                            end_times = data['Observation period'].astype(str).str.split('/').str[1]
                            
                            # Clean up formatting to standard string '2021-01-01 00:01:00'
                            # 将格式清理为标准字符串 '2021-01-01 00:01:00'
                            end_times = end_times.str.replace('T', ' ', regex=False).str.replace('.0', '', regex=False)
                            
                            data = data.drop(columns=['Observation period'])
                            
                            # Round and convert remaining columns to integer types with missing support (Int64)
                            # 四舍五入并将剩余列转换为支持缺失值的整数类型 (Int64)
                            # First round numeric columns / 首先对数值列进行四舍五入
                            num_cols = data.select_dtypes(include=['number']).columns
                            data[num_cols] = data[num_cols].round(0)
                            # Convert radiation columns to Int64 / 将辐射列转换为 Int64
                            for col in num_cols:
                                data[col] = data[col].astype('Int64')
                                
                            # Insert Time explicitly / 显式插入时间列
                            data.insert(0, 'Time', end_times)
                        else:
                            num_cols = data.select_dtypes(include=['number']).columns
                            data[num_cols] = data[num_cols].round(0).astype('Int64')
                        
                        # Ensure csv removes index name 
                        # 确保 csv 移除索引名称
                        data.to_csv(out_file, index=False)
                        break 
                    else:
                        print(f"      Warning: Empty data returned for {stn_id} {yr}")
                        break 
                        
                except Exception as e:
                    print(f"\n      Attempt {attempt+1}/{max_retries} failed for {stn_id} {yr}: {e}")
                    if attempt < max_retries - 1:
                        time.sleep(retry_delay)
                        retry_delay *= 2 
                    else:
                        print(f"      Critical: All {max_retries} attempts failed for {stn_id} {yr}.")
                
                # API Rate limiting 
                # API 频率限制
                time.sleep(2)

if __name__ == "__main__":
    main()

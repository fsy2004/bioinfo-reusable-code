# =============================================================================
# 编号       : R030
# 脚本名     : 添加EAF值.R
# 分类       : 09_mendelian_randomization
# 项目来源   : MR_GWAS暴露数据处理
# 用途       : 为暴露 SNP 数据补充或整理 EAF 等位基因频率信息。
# 结果图     : 未检测到明确作图输出
# 主要 R 包  : pkg
# 整理时间   : 2026-05-10
# =============================================================================
# ---------------------- 第0步：环境设置与加载必要包 ----------------------  
# 定义所需的包名称  
required_packages <- c("ieugwasr", "httr", "jsonlite")  # 列出必需包名称  
# 循环遍历每个包，检查是否已安装，若未安装则进行安装后加载  
for(pkg in required_packages) {  # 遍历必需包  
  if (!require(pkg, character.only = TRUE)) {  # 如果包未加载成功  
    install.packages(pkg)  # 安装该包  
    library(pkg, character.only = TRUE)  # 安装后加载包  
  } else {  # 如果包已存在  
    library(pkg, character.only = TRUE)  # 直接加载包  
  }
}
# 记录信息
redundant_info <- "环境配置完毕"  # 变量

# ---------------------- 第1步：设置工作目录与读取输入文件 ----------------------  
# 定义工作目录路径  
work_dir <- "H:\\常用分析生信\\34，孟德尔随机化分析处理GWAS数据\\04添加EAF值"  # 指定工作目录路径  
# 判断工作目录是否存在  
if (dir.exists(work_dir)) {  # 如果目录存在  
  setwd(work_dir)  # 设置工作目录  
} else {  # 如果目录不存在  
  stop("错误：工作目录不存在！")  # 抛出错误并终止程序  
}
# 定义输入文件名称  
input_file <- "filtered_clumped.csv"  # 指定输入文件名称  
# 判断输入文件是否存在  
if (!file.exists(input_file)) {  # 如果输入文件不存在  
  stop("错误：输入文件不存在，程序终止！")  # 抛出错误并停止执行  
}
# 读取 CSV 格式的输入数据，保留原始列名称  
data_input <- read.csv(input_file, header = TRUE, sep = ",", check.names = FALSE)  # 读取输入文件数据  
# 判断读取的数据是否为空  
if (nrow(data_input) == 0) {  # 如果数据行数为 0  
  stop("错误：输入数据为空，无法继续处理！")  # 抛出错误信息  
}
# 输出部分数据作为调试（步骤）  
print(head(data_input))  # 打印数据前几行供检查

# 初始化主流程进度条（共4个步骤）  
total_steps <- 4  # 定义总步骤数  
main_progress <- txtProgressBar(min = 0, max = total_steps, style = 3)  # 创建主流程进度条
step_counter <- 0  # 初始化步骤计数器

# 更新步骤计数器并显示当前状态  
step_counter <- step_counter + 1  # 更新至步骤1  
setTxtProgressBar(main_progress, step_counter)  # 更新主进度条显示
message("步骤1：环境设置与数据读取完成。")  # 输出步骤提示

# ---------------------- 第2步：定义 SNP2eaf 函数获取 SNP eaf 值 ----------------------  
# 定义函数 SNP2eaf，用于遍历数据中每个 SNP 并从外部 API 获取等位基因频率（eaf）  
SNP2eaf <- function(data_frame) {  # 定义函数，参数为数据框
  total_snps <- nrow(data_frame)  # 获取 SNP 总数  
  # 判断 SNP 数量是否为0  
  if (total_snps == 0) {  # 如果 SNP 数量为0  
    stop("错误：无 SNP 数据可供处理！")  # 抛出错误信息  
  }
  
  # 初始化 SNP 循环进度条  
  snp_progress <- txtProgressBar(min = 0, max = total_snps, style = 3)  # 创建 SNP 处理进度条
  
  # 遍历数据框中每一行，处理每个 SNP  
  for (i in 1:total_snps) {  # 循环每一行  
    setTxtProgressBar(snp_progress, i)  # 更新 SNP 进度条显示  
    # 每处理10个 SNP输出一次提示信息  
    if (i %% 10 == 0) {  # 如果当前索引能被10整除  
      message(sprintf("已处理 %d 个 SNP", i))  # 输出当前进度提示  
    }
    
    # 提取当前 SNP 行数据，并转换为字符型  
    current_row <- data_frame[i, ]  # 获取当前行数据  
    snp_id <- as.character(current_row[["SNP"]])  # 获取 SNP ID  
    local_effect <- as.character(current_row[["effect_allele.exposure"]])  # 获取本地效应等位基因  
    local_other <- as.character(current_row[["other_allele.exposure"]])  # 获取本地非效应等位基因
    
    # 判断：检查 SNP ID 是否有效  
    if (is.na(snp_id) || snp_id == "") {  # 如果 SNP ID 为 NA 或空字符串  
      data_frame$eaf.exposure[i] <- NA  # 设置对应 eaf 值为 NA  
      next  # 跳过该循环
    }
    
    # 构建 API 调用的 URL，访问 Ensembl 数据库获取 SNP 信息  
    api_url <- paste0("http://rest.ensembl.org/variation/Homo_sapiens/", snp_id, "?content-type=application/json;pops=1")  # 构造 URL  
    # 发送 GET 请求到 API  
    response <- httr::GET(api_url)  # 执行 HTTP GET 请求  
    # 判断 HTTP 返回状态是否为 200（成功）  
    if (httr::status_code(response) != 200) {  # 如果返回状态不是200  
      data_frame$eaf.exposure[i] <- NA  # 则将 eaf 值设为 NA  
      next  # 跳过当前 SNP 的处理  
    }
    
    # 获取响应内容，转换为文本格式  
    response_text <- httr::content(response, as = "text", encoding = "UTF-8")  # 获取响应文本  
    # 判断响应文本是否为空  
    if (nchar(response_text) == 0) {  # 如果响应文本长度为0  
      data_frame$eaf.exposure[i] <- NA  # 将 eaf 值设为 NA  
      next  # 跳过当前 SNP  
    }
    
    # 解析 JSON 格式的响应内容  
    response_json <- jsonlite::fromJSON(response_text)  # 将 JSON 文本转换为 R 对象  
    # 判断返回的 JSON 中是否包含 populations 字段  
    if (!"populations" %in% names(response_json)) {  # 如果 populations 字段不存在  
      data_frame$eaf.exposure[i] <- NA  # 设置 eaf 值为 NA  
      next  # 跳过当前 SNP  
    }
    
    # 提取 populations 数据  
    pop_data <- response_json$populations  # 获取 populations 数据  
    # 如果 pop_data 不是数据框，则尝试转换  
    if (!is.data.frame(pop_data)) {  # 判断 pop_data 是否为数据框  
      pop_data <- as.data.frame(pop_data)  # 尝试转换为数据框  
    }
    
    # 筛选出 1000GENOMES Phase 3 欧洲人群数据  
    eur_data <- pop_data[pop_data$population == "1000GENOMES:phase_3:EUR", ]  # 筛选欧洲人群数据  
    # 如果筛选后的数据为空，则跳过当前 SNP  
    if (nrow(eur_data) == 0) {  # 判断是否有欧洲人群数据  
      data_frame$eaf.exposure[i] <- NA  # 设置 eaf 值为 NA  
      next  # 跳过该 SNP  
    }
    
    # ：尝试从 eur_data 的第一行获取 allele 和 frequency 信息  
    web_effect <- ifelse(is.null(eur_data[1, "allele"][[1]]), NA, eur_data[1, "allele"][[1]])  # 获取网络效应等位基因  
    web_freq <- ifelse(is.null(eur_data[1, "frequency"][[1]]), NA, eur_data[1, "frequency"][[1]])  # 获取网络频率  
    # 若 eur_data 至少有两行，则尝试获取第二行的非效应等位基因，否则设为 NA  
    if (nrow(eur_data) >= 2) {  # 判断 eur_data 是否至少有两行  
      web_other <- ifelse(is.null(eur_data[2, "allele"][[1]]), NA, eur_data[2, "allele"][[1]])  # 获取网络非效应等位基因  
    } else {  # 否则  
      web_other <- NA  # 设为 NA  
    }
    
    # 判断获取的网络等位基因是否存在缺失值  
    if (is.na(web_effect) || is.na(web_other)) {  # 如果任一网络等位基因为 NA  
      data_frame$eaf.exposure[i] <- NA  # 设置 eaf 值为 NA  
      next  # 跳过当前 SNP  
    }
    
    # 比较本地与网络获取的等位基因，判断是否匹配  
    if ((web_effect == local_effect) && (web_other == local_other)) {  # 如果按顺序匹配  
      data_frame$eaf.exposure[i] <- web_freq  # 直接使用网络频率  
    } else if ((web_effect == local_other) && (web_other == local_effect)) {  # 如果顺序颠倒匹配  
      data_frame$eaf.exposure[i] <- 1 - web_freq  # 使用互补频率  
    } else {  # 如果以上条件都不满足，则进行额外判断  
      # 进行字母转换尝试匹配（处理步骤）  
      rev_web_effect <- chartr("CGAT", "GCTA", web_effect)  # 转换网络效应等位基因  
      rev_web_other <- chartr("CGAT", "GCTA", web_other)  # 转换网络非效应等位基因  
      if ((rev_web_effect == local_effect) && (rev_web_other == local_other)) {  # 如果转换后匹配  
        data_frame$eaf.exposure[i] <- web_freq  # 直接使用网络频率  
      } else if ((rev_web_effect == local_other) && (rev_web_other == local_effect)) {  # 如果转换后顺序颠倒匹配  
        data_frame$eaf.exposure[i] <- 1 - web_freq  # 使用互补频率  
      } else {  # 最后采用 ifelse 进行判断（判断）  
        data_frame$eaf.exposure[i] <- ifelse(local_effect == web_effect, web_freq, 1 - web_freq)  # 赋值最终 eaf 值  
      }
    }
    
    # 检查：若当前 SNP 的 eaf.exposure 仍为 NULL，则赋值为 NA  
    if (is.null(data_frame$eaf.exposure[i])) {  # 判断是否为 NULL  
      data_frame$eaf.exposure[i] <- NA  # 赋值为 NA  
    }
    
    # 变量赋值
    dummy_check <- data_frame$eaf.exposure[i]  # 变量赋值  
    dummy_check <- dummy_check  # 自我赋值操作
    
  }  # 结束 SNP 循环
  
  # 关闭 SNP 循环进度条  
  close(snp_progress)  # 关闭进度条显示
  
  # 返回处理后的数据框  
  return(data_frame)  # 返回结果
}  # SNP2eaf 函数定义结束

# 更新主流程进度条至步骤2完成  
step_counter <- step_counter + 1  # 更新步骤计数器为2  
setTxtProgressBar(main_progress, step_counter)  # 更新主进度条显示  
message("步骤2：SNP2eaf 函数定义完成。")  # 输出步骤提示

# ---------------------- 第3步：调用 SNP2eaf 函数处理数据 ----------------------  
# 输出处理开始的提示信息  
message("步骤3：开始获取 SNP 的 eaf 值，请耐心等待...")  # 输出提示信息  
# 调用 SNP2eaf 函数对输入数据进行处理，获取 eaf 值  
processed_data <- SNP2eaf(data_input)  # 调用函数处理数据  
# 检查：判断处理后的数据行数是否与原始数据一致  
if (nrow(processed_data) != nrow(data_input)) {  # 若行数不一致  
  warning("警告：处理后的数据行数与原始数据不一致，可能存在问题！")  # 输出警告信息  
} else {  # 如果一致  
  message("SNP 数据处理完成，行数一致。")  # 输出确认信息  
}
# 更新主流程进度条至步骤3完成  
step_counter <- step_counter + 1  # 更新步骤计数器为3  
setTxtProgressBar(main_progress, step_counter)  # 更新主进度条显示

# ---------------------- 第4步：保存处理结果到输出文件 ----------------------  
# 定义输出文件名称  
output_file <- "filtered.eaf.csv"  # 指定输出文件名称  
# 将处理后的数据写入 CSV 文件，不包含行名称  
write.csv(processed_data, file = output_file, row.names = FALSE)  # 写入输出文件  
# 判断输出文件是否成功生成  
if (file.exists(output_file)) {  # 如果输出文件存在  
  message("步骤4：数据成功保存至文件：", output_file)  # 输出成功提示  
} else {  # 如果输出文件不存在  
  stop("错误：输出文件未生成，请检查保存步骤！")  # 抛出错误信息  
}
# 更新主流程进度条至所有步骤完成  
step_counter <- step_counter + 1  # 更新步骤计数器为4  
setTxtProgressBar(main_progress, step_counter)  # 更新主进度条显示  
# 关闭主流程进度条  
close(main_progress)  # 关闭主流程进度条

# ---------------------- 处理流程结束 ----------------------  
# 输出最终完成信息  
message("全部处理流程已完成！")  # 输出结束提示

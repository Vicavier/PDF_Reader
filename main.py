import fitz
import re
import os
from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage,AIMessage
from dotenv import load_dotenv
from tqdm import tqdm
from chapter import CHAPTER_PAGE_RANGE

load_dotenv()

def clean_newlines(s):
    s = re.sub(r'-\n', '\n', s)
    s = re.sub(r'(?<!\.)\n', '', s)
    return s


llm = ChatOpenAI(
    model="deepseek-chat",
    base_url="https://api.deepseek.com/v1",
    api_key=os.getenv("DEEPSEEK_API_KEY")
)
chapter = 2

doc = fitz.open("Coevolutionary Computation and Its Applications.pdf")
print(f"总页数：{doc.page_count}")


with open(CHAPTER_PAGE_RANGE[chapter]["File"],"a",encoding="utf-8") as f:

    for page_num in tqdm(range(CHAPTER_PAGE_RANGE[chapter]["Page"][0],CHAPTER_PAGE_RANGE[chapter]["Page"][1])):
        page = doc.load_page(page_num)
        if page_num == CHAPTER_PAGE_RANGE[chapter]["Page"][0]:
            text_blks = page.get_text("blocks")[:-2]
        else:
            text_blks = page.get_text("blocks")[2:]
        # print(text_blks)
        user_propmt = ""
        sys_prompt = "你是一个演化计算领域的翻译专家，我将给你一段演化计算相关的英文论文段落，请你利用你的专业术语，将英文准确地翻译成中文。注意，区别段落中的标题，公式，等式和数学符号。"
        for text in text_blks:
            user_propmt += f"{clean_newlines(text[4])}"
            
        # 把下一页的第二段接上来
        
        page = doc.load_page(page_num + 1)
        text = page.get_text("blocks")[1]
        user_propmt += f" {text[4]}"

        message = [
            SystemMessage(content=sys_prompt),
            HumanMessage(content=user_propmt)
        ]
        response = llm.invoke(message)
        f.write(response.content)
        f.write("\n")
        f.flush
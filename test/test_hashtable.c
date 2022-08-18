/*
 * @Author: HeYuwei
 * @Date: 2022-08-18 17:50:31
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-08-18 18:12:53
 * @FilePath: \SQP_c\test\test_hashtable.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_TABLE_SIZE 	1024*1024		//最大表长
#define TRUE			1
#define FLASE			0

typedef struct hash_node
{
	struct hash_node *next;	//如果hash(key)相同，依次往后接力
	char *key;				//关键字
	void *value;			//值
	char is_occupyed;	//是否被占用
}Hash_node;

typedef struct hash_table
{
	Hash_node **table;	//哈希表
}Hash_Table;

/*
	功能：		哈希算法，返回哈希值
	key：		char *型，键值
	ken_len:	键值长度
	返回值：	哈希值
*/
static unsigned int JSHash(char* key, unsigned int key_len)
{
	unsigned int hash = 1315423911;
	unsigned int i = 0;

	for(i = 0; i < key_len; key++, i++)
	{
		hash ^= ((hash << 5) + (*key) + (hash >> 2));
	}

	return hash;
}

/*
	功能：		初始化哈希节点
	node：		哈希节点
	返回值：	无
*/
static void init_hs_node(Hash_node *node)
{
	node->next = NULL;
	node->key = NULL;
	node->value = NULL;
	node->is_occupyed = FLASE;
}

/*
	功能：		创建一张哈希表
	参数：		无
	返回值：	成功返回一张建立好的哈希表，失败返回NULL
*/
static Hash_Table *creat_hash_table(void)		//创建一张哈希表
{
	Hash_Table *Hs_table = (Hash_Table *)malloc(sizeof(Hash_Table));	//分配哈希表起始地址
	if (!Hs_table)
	{
		printf("no enough memory\n");
		return NULL;
	}
	
	Hs_table->table = malloc(sizeof(Hash_node) * MAX_TABLE_SIZE);	//为哈希表的所有存储节点分配内存
	if (!Hs_table->table)
	{
		printf("no enough memory for table\n");
		free(Hs_table);
		Hs_table = NULL;
		return NULL;
	}
	
	memset(Hs_table->table, 0, sizeof(Hash_node) * MAX_TABLE_SIZE);	//所有节点初始化为0
	
	return Hs_table;		//返回哈希表首地址
}

/*
	功能：		向哈希表中增加一个节点
	Hs_table：	哈希表，不能为空
	key：		键值，不能为空
	key_len:	键值长度
	value:		值
	返回值：	0 成功   -1 失败
*/
int add_node2HashTable(Hash_Table *Hs_table, char *key, unsigned int key_len, void *value)
{
	if(!Hs_table || !key )
	{
		printf("something is NULL\n");
		return -1;
	}

	unsigned int i = JSHash(key, key_len) % MAX_TABLE_SIZE;		//通过jhash取其键值对应的哈希表下标

	Hash_node *p = Hs_table->table[i];
	Hash_node *pri = p;
	
	while(p)	//该点已经有哈希节点了，需要走到这条链的最后
	{
		if ( strncmp(key, p->key, key_len) == 0 )	//该键值已经存在，更新value
		{
			if(p->is_occupyed)
			{
				p->value = value;
				p->is_occupyed = 1;		//占用标记置1
				break;
			}
		}
		pri = p;			//pri始终指向p的上一个位置，用于保留最后一个table[i]的最后一个哈希节点
		p = p->next;
	}
	
	if(!p)	//走到最后或者该点一直没有被占用
	{
		Hash_node *tmp = (Hash_node *)malloc(sizeof(Hash_node));
		if( !tmp )
		{
			printf("no enough memory\n");
			return -1;
		}
		init_hs_node(tmp);
		char *tmp_key = (char *)malloc(key_len+1);
		if(!tmp_key)
		{
			free(tmp);
			tmp = NULL;
			return -1;
		}
		strncpy(tmp_key, key, key_len);
		tmp->key = tmp_key;
		tmp->value = value;
		tmp->is_occupyed = TRUE;	//更新占用标记
		
		if(pri == NULL)		//该点没有被占用过，直接指
		{
			Hs_table->table[i] = tmp;
		}
		else		//该点被占用过，已经走到这条链最后
		{
			pri->next = tmp;
		}
	}
	
	return 0;
}

/*
	功能：		从哈希表中获取数据
	Hs_table：	哈希表，不能为空
	key：		键值，不能为空
	key_len:	键值长度
	返回值：	存储内容，未找到则为NULL
*/
void *get_value_from_hstable(Hash_Table *Hs_table, char *key, unsigned int key_len)
{
	if( !Hs_table || !key)
	{
		printf("something is NULL\n");
		return NULL;
	}
	
	int i = JSHash(key,key_len) % MAX_TABLE_SIZE;
	Hash_node *tmp = Hs_table->table[i];
	
	while(tmp)
	{
		if(strncmp(tmp->key, key, key_len) == 0)
		{
			return tmp->value;
		}
		tmp = tmp->next;
	}
	
	return NULL;
}

/*
	功能：		删除哈希表
	Hs_Table:	哈希表名称
	返回值：	无
*/
void hash_table_delete(Hash_Table *Hs_Table)
{
    if (Hs_Table)
	{
        if (Hs_Table->table)
		{
            int i = 0;
            for (i = 0; i<MAX_TABLE_SIZE; i++)	//遍历每一个table 
			{
                Hash_node *p = Hs_Table->table[i];
                Hash_node *q = NULL;
                while (p)		//该点存在存储内容
				{
                    q = p->next;
                    p->is_occupyed = 0;	//占用标记清0
                    p = q;
                }
            }
            free(Hs_Table->table);	//释放表存储指针的内存占用
            Hs_Table->table = NULL;
        }
        free(Hs_Table);		//释放表指针
		Hs_Table = NULL;
    }
	return ;
}

typedef struct commodity
{
	char name[32];	//商品名称
	float price;	//价格
}Com;

void printf_com_info(Com *com)
{
	printf("name=%s\tprice=%.1f\n", com->name, com->price);
	return ;
}

int main(int argc, char **argv)
{
	Hash_Table *Hs_table = creat_hash_table();
	if(!Hs_table)
	{
		printf("creat hash table fail\n");
		return -1;
	}
	char name[32] = {0};
	int i = 0;
	for (i = 0; i < 100; i++)
	{
		Com *tmp_com = (Com *)malloc(sizeof(Com));
		//先不对内存分配进行判断
		sprintf(tmp_com->name, "com%d", i);
		tmp_com->price = rand() % 1000;
		
		add_node2HashTable(Hs_table, tmp_com->name, strlen(tmp_com->name), tmp_com);
	}

	for(i = 0; i < 100; i++)
	{
		memset(name, 0, sizeof(name));
		sprintf(name, "com%d", i);
		Com *get_com = get_value_from_hstable(Hs_table, name, strlen(name));
		if(get_com)
		{
			printf_com_info(get_com);
		}
	}

	hash_table_delete(Hs_table);

	return 0;
}

